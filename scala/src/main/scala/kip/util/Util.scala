package kip.util

import scala.util.matching.Regex
import java.io.{BufferedWriter, FileWriter}



object Util {
  case class Vec3(x: Double, y: Double, z: Double) {
    def distance2(that: Vec3) = {
      sqr(x-that.x) + sqr(y-that.y) + sqr(z-that.z)
    }
  }
  
  final def sqr(x: Double) = x*x
  
  // def refToOption[A <: AnyRef](a:A) =
  //  if (a == null) None else Some(a)  

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
}
