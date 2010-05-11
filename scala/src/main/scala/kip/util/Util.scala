package kip.util

import scala.util.matching.Regex
import java.io.{LineNumberReader, BufferedReader, BufferedWriter, FileWriter, FileReader}



object Util {
  case class Vec3(x: Double, y: Double, z: Double) {
    def distance2(that: Vec3) = {
      sqr(x-that.x) + sqr(y-that.y) + sqr(z-that.z)
    }
  }
  
  final def sqr(x: Double) = x*x
  
  final def printTime[A](fn: => A) = {
    val begin = System.currentTimeMillis
    print ("Begin timing... ")
    val ret = fn
    val end = System.currentTimeMillis
    println("done. Execution time: " + (end-begin)/1000. + "s")
    ret
  }

  def average(vs: Seq[Double]) =
    if (vs.size == 0) Double.NaN else vs.reduceLeft(_+_) / vs.length
  
  def trace[A](str: String, v: A) = {
    println(str + v)
    v
  }

  def formatDataInColumns(kvs: (String, Array[Double])*) = {
    val sb = new StringBuffer()

    val (descs, vals) = kvs.unzip
    sb.append("# " + descs.mkString(" ") + "\n")
    for (i <- 0 until vals(0).size) {
      sb.append(vals.map{_(i)}.mkString(" ")+"\n")
    }
    sb.toString()
  }

  def writeStringToFile(s: String, fn: String) {
    val writer = new BufferedWriter(new FileWriter(fn))
    writer.write(s);
    writer.close();
  }
  
  def readDataFromFile(fn: String): (Array[Array[Double]], String) = {
    import scala.collection.mutable.ArrayBuffer
    var data: Array[ArrayBuffer[Double]] = null
    var desc: String = ""
    
    val reader = new LineNumberReader(new FileReader(fn))
    var line = reader.readLine()
    while (line != null) {
      line = line.trim
      
      if (line.size == 0) {
        // skip
      }
      else if (line.startsWith("#")) {
        desc += line
      }
      else {
        val items = line.split("\\s+").map(_.toDouble)
        if (data == null) {
          data = Array.fill(items.size)(new ArrayBuffer[Double]())
        }
        if (data.size != items.size) {
          throw new Error("Column size mismatch, "+data.size+"!="+items.size+", in '"+fn+"':"+reader.getLineNumber())
        }
        for (i <- items.indices) {
          data(i).append(items(i))
        }
      }
      line = reader.readLine()
    }
    
    (data.map(_.toArray), desc) 
  }

}
