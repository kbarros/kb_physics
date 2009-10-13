package kip

import scala.io.Source
import scala.collection.mutable.ArrayBuffer
import scikit._
import Util._



case class Vec3(x: Double, y: Double, z: Double) {
    final def distance2(that: Vec3) =
	sqr(x-that.x) + sqr(y-that.y) + sqr(z-that.z) 
}


class Snapshot(val time: Double, val natoms: Int) {	
    var lo: Vec3 = null
    var hi: Vec3 = null
    
    var id : Array[Double] = null
    var typ : Array[Double] = null
    var x : Array[Double] = null
    var y : Array[Double] = null
    var z : Array[Double] = null
    var ix : Array[Double] = null
    var iy : Array[Double] = null
    var iz : Array[Double] = null
    var vx : Array[Double] = null
    var vy : Array[Double] = null
    var vz : Array[Double] = null
    
    var temperature = Double.NaN
    var energy = Double.NaN
    var pressure = Double.NaN
    
    def getPoint(i: Int) = 
	Vec3(x(i), y(i), z(i))
}


object LammpsParser {	
    
    // Similar to String.split() but marginally faster
    /*
    def splitWhitespace(line: String) = {
	val tokens = new ArrayBuffer[String]
	val chars = line.toCharArray()
	
	def isWhitespace(c: Char) = c == ' '
	
	var offset = 0
	while (offset < chars.length) {
	    while (offset < chars.length && isWhitespace(chars(offset)))
		offset += 1
	    var count = 0
	    while (offset+count < chars.length && !isWhitespace(chars(offset+count)))
		count += 1
	    tokens.append(new String(chars, offset, count))
	    offset += count
	}
	tokens
    }
    */
    
    def readSnapshot(lines: Iterator[String]) = {
	def matchLine(line: String, str: String) {
	    if (line != str)
		throw new Exception("Failure to parse '"+str+"'")
	}
	
	matchLine(lines.next, "ITEM: TIMESTEP\n")
	val time = lines.next.toDouble
	
	matchLine(lines.next, "ITEM: NUMBER OF ATOMS\n")		
	val natoms = lines.next.trim.toInt

	matchLine(lines.next, "ITEM: BOX BOUNDS\n")
	val Array(xlo:Double,xhi) = lines.next.split("\\s").map(_.toDouble)
	val Array(ylo:Double,yhi) = lines.next.split("\\s").map(_.toDouble)
	val Array(zlo:Double,zhi) = lines.next.split("\\s").map(_.toDouble)
	
	// read list of atoms
	val descLine = lines.next
	if (!descLine.startsWith("ITEM: ATOMS "))
	    throw new Exception("Failure to parse 'ITEMS: ATOMS'")
	val descs = descLine.split("\\s").drop(2)
	val cols = new Array[Array[Double]](descs.length, natoms)
	for (atom <- 0 until natoms) {
	    val vs = lines.next.split("\\s")
//	    val vs = splitWhitespace(lines.next)
	    for (value <- 0 until vs.length) {
		cols(value)(atom) = vs(value).toDouble 
	    }
	}
	
	// build snapshot object
	val ss = new Snapshot(time, natoms)
	ss.lo = Vec3(xlo, ylo, zlo)
	ss.hi = Vec3(xhi, yhi, zhi)
	for ((d,c) <- descs zip cols) {
	    d match {
		case "id" => ss.id = c // atom index, using Fortran indexing convention
		case "type" => ss.typ = c
		case "x" => ss.x = c
		case "y" => ss.y = c
		case "z" => ss.z = c
		case "ix" => ss.ix = c
		case "iy" => ss.iy = c
		case "iz" => ss.iz = c
		case "vx" => ss.vx = c
		case "vy" => ss.vy = c
		case "vz" => ss.vz = c
	    }
	}
	
	// wrap particle positions according to image indices (ix, iy, iz)
	if (ss.ix != null) {
	    for (i <- 0 until ss.x.length) {
		ss.x(i) += ss.ix(i)*(xhi-xlo)
		ss.y(i) += ss.iy(i)*(yhi-ylo)
		ss.z(i) += ss.iz(i)*(zhi-zlo)
	    }
	}
	
	ss
    }
    
    def readLammpsDumpPartial(fname: String, maxSnapshots: Int) = {
	val lines = Source.fromFile(fname).getLines
	val snaps = new ArrayBuffer[Snapshot]
	while (lines.hasNext && snaps.length < maxSnapshots) {
	    snaps.append(readSnapshot(lines))
	    if (snaps.length % 1000 == 0)
		System.err.println("Reading snapshot "+snaps.length)
	}
	snaps
    }
    
    def readLammpsDump(fname: String) = {
	readLammpsDumpPartial(fname, Integer.MAX_VALUE)
    }
/*
    def readLammpsThermo(fname: String, snaps: Seq[Snapshot]) {
	val lines = Source.fromFile(fname).getLines
	val snaps = new ArrayBuffer[Snapshot]
	
	var desc = None
	def nextThermoLine = {
	    // parse until next descriptor
	    if (desc == None) {
		match (lines.find(_.startsWith("Step "))) {
		    case None => return None
		    case Some(d) => desc = d
		}
	    }
	    
	    if (!lines.hasNext)
		return None
	}
	
	for (s <- snaps) {
	    parse
	}
	var idx = 0
	while (idx < runIndex) {
	    if ()
	}
	while (lines.hasNext) {
	    snaps.append(readSnapshot(lines))
	    if (snaps.length % 1000 == 0)
		System.err.println("Reading snapshot "+snaps.length)
	}
	snaps
    }
    */
}

