package kip.util

import java.io.{BufferedReader, File, FileInputStream, FileReader, InputStreamReader}
import java.util.zip.GZIPInputStream
import scala.io.Source
import scala.collection.mutable.ArrayBuffer

import kip.math.Vec3
import kip.math.Math.sqr


class Thermo {
  var time = Integer.MIN_VALUE
  var temperature = Double.NaN
  var energy = Double.NaN
  var pressure = Double.NaN
}


class Snapshot(val time: Int, val natoms: Int) {
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
  var q  : Array[Double] = null
  
  var thermo: Thermo = null

  def getPoint(i: Int) = Vec3(x(i), y(i), z(i))
  
  def distance2(i: Int, j: Int) = {
    val lx = hi.x - lo.x
    val ly = hi.y - lo.y
    val lz = hi.z - lo.z
    var dx = math.abs(x(i) - x(j))
    var dy = math.abs(y(i) - y(j))
    var dz = math.abs(z(i) - z(j))
    while (dx > lx/2) dx -= lx
    while (dy > ly/2) dy -= ly
    while (dz > lz/2) dz -= lz
    dx*dx + dy*dy + dz*dz
  }

  def maxDistance2 = {
    val lx = hi.x - lo.x
    val ly = hi.y - lo.y
    val lz = hi.z - lo.z
    sqr(lx/2) + sqr(ly/2) + sqr(lz/2)
  }
  
  def volume = {
    val lx = hi.x - lo.x
    val ly = hi.y - lo.y
    val lz = hi.z - lo.z
    lx * ly * lz
  }
}

object LammpsParser {  
  def readLammpsThermo(fname: String) = {
    val lines = Source.fromFile(new File(fname)).getLines().buffered
    
    val dataLines = new ArrayBuffer[Array[Double]]()
    var desc = None: Option[Array[String]]
    
    while (lines.hasNext) {
      var line = lines.next
      if (line.startsWith("Step ")) {
        desc = desc.orElse(Some(line.split("\\s").filter(_.length > 0)))
        
        while (lines.hasNext && !lines.head.startsWith("Loop time")) {
          line = lines.next
          dataLines.append(line.split("\\s").filter(_.length > 0).map(_.toDouble))
        }
      }
    }
    
    dataLines map { vs =>
      var thermo = new Thermo
      for ((d,v) <- desc.get zip vs) {
        d match {
          case "Step" => thermo.time = v.toInt
          case "Temp" => thermo.temperature = v
          case "TotEng" => thermo.energy = v
          case "Press" => thermo.pressure = v
          case _ => ()
        }
      }
      thermo
    }
  }
  
  def readSnapshot(lines: Iterator[String]) = {
    def matchLine(line: String, str: String) {
      if (line != str)
        throw new Exception("Failure to parse '"+str+"'; found '"+line+"'")
    }

    matchLine(lines.next, "ITEM: TIMESTEP")
    val time = lines.next.toInt

    matchLine(lines.next, "ITEM: NUMBER OF ATOMS")
    val natoms = lines.next.trim.toInt

    matchLine(lines.next, "ITEM: BOX BOUNDS")
    val Array(xlo, xhi) = lines.next.split("\\s").map(_.toDouble)
    val Array(ylo, yhi) = lines.next.split("\\s").map(_.toDouble)
    val Array(zlo, zhi) = lines.next.split("\\s").map(_.toDouble)

    // read list of atoms
    val descLine = lines.next
    if (!descLine.startsWith("ITEM: ATOMS "))
      throw new Exception("Failure to parse 'ITEMS: ATOMS'")
    val descs = descLine.split("\\s").drop(2)
    val cols = Array.ofDim[Double](descs.length, natoms)
    val vs = Array.ofDim[Double](descs.length)
    for (atom <- 0 until natoms) {
      val cnt = kip.util.Parsing.stringSplitDouble(vs, lines.next)
      assert(cnt == descs.length)
      for (value <- 0 until descs.length) {
        cols(value)(atom) = vs(value)
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
        case "q"  => ss.q = c
        case "fx" => ()
        case "fy" => ()
        case "fz" => ()
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

  def readLammpsDump(
    fname: String,
    process: Snapshot => Option[Snapshot] = Some(_),
    maxSnapshots: Int = Integer.MAX_VALUE
  ): Seq[Snapshot] = {
    val isGZipped = fname.endsWith("gz") || fname.endsWith("gzip")
    val fis = new FileInputStream(fname)
    val gzis = if (isGZipped) new GZIPInputStream(fis) else fis
    val br = new BufferedReader(new InputStreamReader(gzis))
    val lines = Iterator.continually {
      try { br.readLine() } catch {
        case e: java.io.EOFException => {
          System.err.println("Unexpected end of file '"+fname+"'")
          null
        }
      }
    } takeWhile (_ != null)
    
    val snaps = new ArrayBuffer[Snapshot]()

    try {
      while (lines.hasNext && snaps.length < maxSnapshots) {
        //if (snaps.length % 100 == 0)
        //  System.err.println("Reading snapshot "+snaps.length)
        process(readSnapshot(lines)) foreach { snaps.append(_) }
      }
    } catch {
      case e: Exception => {
        System.err.println(e.toString)
        System.err.println("Failed to parse snapshot "+(snaps.size+1))
      }
    }
    
    snaps
  }
  
  def weaveThermoData(snaps: Seq[Snapshot], thermos: Seq[Thermo]) {
    val thermoMap = new scala.collection.mutable.HashMap[Int,Thermo]
    for (th <- thermos) {
      thermoMap += (th.time -> th)
    }
    for (sn <- snaps) {
      sn.thermo = thermoMap.get(sn.time).getOrElse(new Thermo)
    }
  }
  
  def main(args: Array[String]) {
    System.out.println("hello")
    readLammpsThermo("/Users/kbarros/Desktop/log.lammps") foreach {th => println(th.time)}
  }
}

