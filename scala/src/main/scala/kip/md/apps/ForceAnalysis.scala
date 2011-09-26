package kip.md.apps

object ForceAnalysis extends App {
  val lines = scala.io.Source.fromFile(args(0)).getLines
  
  val entries: Array[(Double, Double)] = lines.toArray map { s =>
    val Array(t, f) = s.split("\\s")
    (t.toDouble, f.toDouble)
  }
  
  // calculate
  val minTime = 400.0
  val scale = 2.0 // convert force to stress
  
  val entriesCut = entries.dropWhile(_._1 < minTime)
  val ts = entriesCut.map(_._1)
  val fs = entriesCut.map(_._2 * scale)
  
  import kip.util.Statistics._
  println("average = %g".format(mean(fs)))
  println("stdev   = %g".format(stddev(fs)))
  
  val fk = kip.math.Math.fftReal(fs).map(_.abs)
  val ks = fk.indices.map(_.toDouble).toArray
  
  val plot = new scikit.graphics.dim2.Plot("Force")
  scikit.util.Utilities.frame(plot.getComponent(), plot.getTitle())
  plot.registerLines("Data", new scikit.dataset.PointSet(ks drop 1, fk drop 1), java.awt.Color.BLACK)
}
