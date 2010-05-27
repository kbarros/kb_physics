package kip.util

import scala.util.matching.Regex
import java.io.{LineNumberReader, BufferedReader, BufferedWriter, FileWriter, FileReader}



object Util {
  case class Vec3(x: Double, y: Double, z: Double) {
    def distance2(that: Vec3) = {
      sqr(x-that.x) + sqr(y-that.y) + sqr(z-that.z)
    }
  }
  
  def sqr(x: Double) = x*x
  
  def average(vs: Seq[Double]) =
    if (vs.size == 0) Double.NaN else vs.sum / vs.length

  def time[A](f: => A, s: String): A = {
    print("Timing '"+s+"'...")
    val t1 = System.currentTimeMillis
    val ret = f
    val t2 = System.currentTimeMillis
    println(" done. Elapsed time "+(t2-t1)/1000.+"s");
    ret
  }
  
  def tr[A](v: A, str: String) = {
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
