package kip.projects.dlc

import math._

object IntegratedNeighbors extends App {
  
  
  def integratePairCorrelation(r: Array[Double], gr: Array[Double]): Array[Double] = {
    val ret = new Array[Double](r.size)
  
    val L = 72.3279822728
    val rho = 100.0 / (L*L*L)

    var acc = 0.0

    for (i <- 0 until r.size-1) {
      val dr = r(i+1) - r(i)
      acc += 4 * Pi * r(i) * r(i) * dr * gr(i) * rho
      ret(i) = acc
    }
    ret
  }
  
  val data = args.map(fn => kip.util.Util.readCsvFile(fn)._1)
  
  val r = data(0)(0)
  val grs = data.map(d => integratePairCorrelation(r, d(1)))
  
  
  println("# r <neighbors>* ")
  for (i <- 0 until r.size-1) {
    val row = r(i) +: grs.map(gr => gr(i))
    println(row.mkString(" "))
  }
  println()
}
