package kip

import java.lang.Math._
import java.io.{BufferedWriter, FileWriter}
import scala.io.Source
import scala.collection.mutable.ArrayBuffer
import scikit._


object Util {
    final def printTime[A](fn: => A) = {
	val begin = System.currentTimeMillis
	print ("Begin timing... ")
	val ret = fn
	val end = System.currentTimeMillis
	println("done. Execution time: " + (end-begin)/1000. + "s")
	ret
    }
    
    final def sqr(x: Double) =
	x*x
    
    def average(vs: Seq[Double]) =
	vs.reduceLeft(_+_) / vs.length
    
    def trace[A](str: String, v: A) = {
	println(str + v)
	v
    }
    
    def fft1d_continuous(a: Array[Double], L:Double) = {
	import scikit.numerics.fft.FFT1D
	import scikit.dataset.Accumulator

	val accumulator = new Accumulator()
	
	val fft = new FFT1D(a.length)
	fft.setLength(L)
	fft.transform(a, new FFT1D.MapFn() {
	    def apply(k:Double, re:Double, im:Double) {
		accumulator.accum(k, re*re+im*im)
	    }
	})
	accumulator
    }
    
    def formatDataInColumns(desc: String, vs: Seq[Array[Double]]) = {
	val sb = new StringBuffer()
	sb.append(desc + "\n")
	for (i <- 0 until vs(0).size) {
	    vs.foreach {v => sb.append(v(i)+" ")}
	    sb.append("\n")
	}
	sb.toString()
    }

    def writeStringToFile(s: String, fn: String) {
	val writer = new BufferedWriter(new FileWriter(fn))
    	writer.write(s);
    	writer.close();
    }
    
    def plot(a: Array[Double]) {
	scikit.util.Commands.plot(a)
    }
}

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
    
    
    def getPoint(i: Int) = 
	Vec3(x(i), y(i), z(i))
    
    def centerOfMass(elems: Seq[Int]) = {
	var xAcc = 0.
	var yAcc = 0.
	var zAcc = 0.
	for (i <- elems) {
	    xAcc += x(i)
	    yAcc += y(i)
	    zAcc += z(i)
	}
	val n = elems.length
	Vec3(xAcc/n, yAcc/n, zAcc/n)
    }
    
    def radiusOfGyrationSquared(elems: Seq[Int]) = {
	val cm = centerOfMass(elems)
	
	var x2Acc = 0.
	var y2Acc = 0.
	var z2Acc = 0.
	for (i <- elems) {
	    x2Acc += sqr(x(i)-cm.x) 
	    y2Acc += sqr(y(i)-cm.y) 
	    z2Acc += sqr(z(i)-cm.z)
	}
	val n = elems.length
	(x2Acc/n, y2Acc/n, z2Acc/n)		
    }
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
    
    def matchLine(line: String, str: String) {
	if (line != str)
	    throw new Exception("Failure to parse '"+str+"'")
    }
    
    def readSnapshot(lines: Iterator[String]) = {
	matchLine(lines.next, "ITEM: TIMESTEP\n")
	val time = lines.next.toDouble
	
	matchLine(lines.next, "ITEM: NUMBER OF ATOMS\n")		
	val natoms = lines.next.trim.toInt

	matchLine(lines.next, "ITEM: BOX BOUNDS\n")
	val Array(xlo:Double,xhi) = lines.next.split("\\s").map(_.toDouble)
	val Array(ylo:Double,yhi) = lines.next.split("\\s").map(_.toDouble)
	val Array(zlo:Double,zhi) = lines.next.split("\\s").map(_.toDouble)
	
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
	
	val ss = new Snapshot(time, natoms)
	ss.lo = Vec3(xlo, ylo, zlo)
	ss.hi = Vec3(xhi, yhi, zhi)
	
	for ((d,c) <- descs zip cols) {
	    d match {
		case "id" => () // ss.id = c // redundant information, using Fortran indexing convention
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

	ss
    }
    
    def readLammpsDumpPartial(fname: String, maxSnapshots: Int) = {
	val src = Source.fromFile(fname)
	val lines = src.getLines
	
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
    
    def writeVelocities(snaps: Seq[Snapshot]) {
	val sb = new StringBuffer
	sb.append("# time ")
	for (i <- 0 until snaps(0).natoms) {
	    sb.append("atom"+i+" ")
	}
	sb.append("\n")
	
	for (ss <- snaps) {
	    sb.append (ss.time + " ")
	    for (vx <- ss.vx) {
		sb.append(vx + " ")
	    }
	    sb.append("\n")
	}
	
	writeStringToFile(sb.toString, "/Users/kbarros/Desktop/blah.dat")
    }
}

