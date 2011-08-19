package kip.projects.quantum

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
  
  def plotLines(plot: Plot, data: (Array[R], Array[R])) {
    val pts = new scikit.dataset.PointSet(data._1, data._2)
    plot.registerBars("Coarse", pts, java.awt.Color.RED)
  }
  
  def eigenvaluesExact(H: PackedSparse[S]): Array[R] = {
    val (v, w) = H.toDense.eig
    v.map(_.re).toArray.sorted
  }
  
  def eigenvaluesApprox(H: PackedSparse[S], order: Int): (Array[R], Array[R]) = {
    val kpm = new KPM(H, order)
    kpm.eigenvaluesApprox(kpm.jacksonKernel)
  }
}

// Kernel polynomial method

class KPM(H: PackedSparse[S], order: Int) {
  val n = H.numRows
  
  // calculate moments
  def moments: Array[R] = {
    val ret = Array.fill(order)(0d)
    
    val v0 = dense(n, 1)
    val v1 = dense(n, 1)
    val v2 = dense(n, 1)
    
    // basis vectors
    for (b <- 0 until n) {
      v0.clear(); v0(b) = 1  // alpha_0 = T_0[H] |b> =   |b>
      v1 :=* (H, v0)         // alpha_1 = T_1[H] |b> = H |b>
      
      ret(0) += v0(b).re
      ret(1) += v1(b).re
      
      for (n <- 2 until order) {
        // v0 = alpha_{m-2}
        // v1 = alpha_{m-1}
        // v2 = alpha_m = 2 H v1 - v0
        v2 := v0
        v2.gemm(2, H, v1, -1)
        
        ret(n) += v2(b).re
        v0 := v1
        v1 := v2
      }
    }
    
    ret
  }
  
  lazy val jacksonKernel: Array[R] = {
    import math._
    val Np = order+1
    Array.tabulate(order) { n => 
      (1/Np)*((Np-n)*cos(Pi*n/Np) + sin(Pi*n/Np)/tan(Pi/Np))
    }
  }
  
  def reconstruct(mu: Array[R], g: Array[R], x: R): R = {
    import math._
    var acc = g(0)*mu(0)
    var t0 = x // T_{-1}
    var t1 = 1d // T_{0}
    for (n <- 1 to order-1) {
      // t0 = T_{n-2}
      // t1 = T_{n-1}
      // t2 = T_n
      val t2 = 2*x*t1 - t0
      acc += 2*g(n)*mu(n)*t2
      
      t0 = t1
      t1 = t2
    }
    (1 / Pi * sqrt(1 - x*x)) * acc
  }
  
  def eigenvaluesApprox(kernel: Array[R]): (Array[R], Array[R]) = {
    val mu = moments
    val xs = Array.tabulate[R](order) { i =>
      i / order.toDouble
    }
    val ret = (xs, xs.map(reconstruct(mu, kernel, _)))
    ret
  }
}
