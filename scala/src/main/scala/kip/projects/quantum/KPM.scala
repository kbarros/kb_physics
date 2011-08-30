package kip.projects.quantum

import kip.util.Util.time
import kip.math.Math._
import kip.enrich._
import smatrix._

import ScalarTyp._
import ctor._

object KPM {
  import scikit.graphics.dim2._
  
  def mkPlot(): Plot = {
    val plot = new Plot("Histogram")
    scikit.util.Utilities.frame(plot.getComponent(), plot.getTitle())
    plot
  }
  
  def plotHistogram(plot: Plot, a: Array[R]) {
    val hist = new scikit.dataset.Histogram(0.01)
    a.foreach { hist.accum(_) }
    plot.registerBars("Density", hist, java.awt.Color.BLACK)
  }
  
  def plotLines(plot: Plot, data: (Array[R], Array[R]), name: String, color: java.awt.Color) {
    val pts = new scikit.dataset.PointSet(data._1, data._2)
    plot.registerLines(name, pts, color)
  }
  
  def chebyshev(m: Int, x: Double): Double = {
    if (m == 0) {
      1d
    }
    else {
      var t0 = x // T_{-1}
      var t1 = 1d // T_{0}
      for (i <- 1 to m) {
        val t2 = 2*x*t1 - t0
        t0 = t1 // T_{i-1}
        t1 = t2 // T_i
      }
      t1
    }
  }
  
  // \int dx x^{moment} y(x) 
  def integrate(xs: Array[Double], ys: Array[Double], moment: Int): Array[Double] = {
    require(xs.size >= 2 && xs.size == ys.size)
    require((1 until xs.size).forall(i => xs(i-1) < xs(i)))
    xs.indices.toArray.foldMapLeft(0d) { (acc, i) =>
      val x0 = if (i > 0)         xs(i-1) else xs(i)
      val x2 = if (i < xs.size-1) xs(i+1) else xs(i)
      val x = 0.5*(x0 + x2)
      val dx = 0.5 * (x2 - x0)
      acc + dx*math.pow(x, moment)*ys(i)
    }
  }
  
  def integrateDeltas(xs: Array[Double], deltas: Array[Double], moment: Int): Array[Double] = {
    require(xs.size >= 2)
    require((1 until xs.size).forall(i => xs(i-1) < xs(i)))
    require((1 until deltas.size).forall(i => deltas(i-1) <= deltas(i)))
    require(deltas.head >= xs.head && deltas.last <= xs.last)
    val ret = new Array[Double](xs.length)
    var j = 0
    var acc = 0d
    for (i <- xs.indices) {
      val binEnd = if (i < xs.size-1) avg(xs(i), xs(i+1)) else xs(i)
      while (j < deltas.size && deltas(j) <= binEnd) {
        acc += math.pow(deltas(j), moment)
        j += 1
      }
      ret(i) = acc
    }
    ret
  }
}

// Kernel polynomial method

class KPM(H: PackedSparse[S], order: Int, nrand: Int, seed: Int = 0) {
  val rand = new util.Random(seed)
  val n = H.numRows
  
  
  def momentsStochastic(): Array[R] = {
    val ret = Array.fill(order)(0d)
    
    val r  = dense(n, nrand)
    val t0 = dense(n, nrand)
    val t1 = dense(n, nrand)
    val t2 = dense(n, nrand)
    
    // r random vectors with uniformly distributed values (1, i, -1, -i) 
    r.fill(Seq[S#A](1, I, -1, -I).apply(rand.nextInt(4)))
    
    t0 := r
    t1 :=* (H, t0)         // T_1[H] |r> = H |r>
    
    // Set first two moments are set to their exact values
    ret(0) = n
    ret(1) = H.trace.re
    // ret(0) += (r dagDot t0).re
    // ret(1) += (r dagDot t1).re
    
    for (m <- 2 until order) {
      // t0 = T_{m-2}[H] r
      // t1 = T_{m-1}[H] r
      // t2 = T_m[H] r = 2 H t1 - t0
      t2 := t0
      t2.gemm(2, H, t1, -1)

      ret(m) = (r dagDot t2).re / nrand
      t0 := t1
      t1 := t2
    }
    
    ret
  }
  
  // slow calculation of moments
  def momentsSlow(): Array[R] = {
    val ret = Array.fill(order)(0d)
    
    val Hd = H // H.toDense
    
    val t0 = eye(n).toDense
    val t1 = H.toDense
    val t2 = dense(n, n)
    
    ret(0) = t0.trace.re
    ret(1) = t1.trace.re
    
    for (m <- 2 until order) {
      // t0 = T_{m-2}(H)
      // t1 = T_{m-1}(H)
      // t2 = T_m(H) = 2 H t1 - t0
      t2 := t0
      t2.gemm(2, Hd, t1, -1)
      
      ret(m) = t2.trace.re
      t0 := t1
      t1 := t2
    }
    
    ret
  }

  lazy val jacksonKernel: Array[R] = {
    import math._
    val Np = order+1d
    Array.tabulate(order) { n => 
      (1/Np)*((Np-n)*cos(Pi*n/Np) + sin(Pi*n/Np)/tan(Pi/Np))
    }
  }
  
  def reconstruct(mu: Array[R], g: Array[R], x: R): R = {
    import math._
    var acc = g(0)*mu(0)
    var t0 = x // T_{-1}
    var t1 = 1d // T_{0}
    
    for (m <- 1 to order-1) {
      // t0 = T_{n-2}
      // t1 = T_{n-1}
      // t2 = T_n
      val t2 = 2*x*t1 - t0
      acc += 2*g(m)*mu(m)*t2
      
      t0 = t1
      t1 = t2
    }
    (1 / (Pi * sqrt(1 - x*x))) * acc
  }
  
  val range: Array[R] = {
    val nx = 5*order 
    Array.tabulate[R](nx) { i =>
      2d * i / nx.toDouble - 1d + (0.5 / nx)
    }
  }
  
  def eigenvaluesExact(): Array[R] = {
    val (v, w) = time("Exact diagonalization of N=%d matrix".format(H.numRows))(H.toDense.eig)
    v.map(_.re).toArray.sorted
  }
  
  def eigenvaluesApprox(kernel: Array[R]): Array[R] = {
//    val mu = time("Calculating %d moments of N=%d matrix".format(order, H.numRows))(moments())
//    val mu = time("Calculating %d moments (slow)".format(order))(momentsSlow())
    
    val mu = time("Calculating %d moments (stoch) of N=%d matrix".format(order, H.numRows))(momentsStochastic())
      
    time("Reconstructing density")(range.map(reconstruct(mu, kernel, _)))
  }
}
