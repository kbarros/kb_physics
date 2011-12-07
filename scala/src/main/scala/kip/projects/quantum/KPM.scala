package kip.projects.quantum

import kip.util.Util.{time, notime}
import kip.math.Math._
import kip.enrich._
import smatrix._
import ctor._


object KPM {
  import scikit.graphics.dim2._

  def chebyshevFillArray(x: Double, ret: Array[Double]) {
    if (ret.size > 0)
      ret(0) = 1
    if (ret.size > 1)
      ret(1) = x
    for (m <- 2 until ret.size) {
      ret(m) = 2*x*ret(m-1) - ret(m-2)
    }
  }
  
  def chebyshevArray(x: Double, m: Int): Array[Double] = {
    val ret = Array.fill[Double](m)(0d)
    chebyshevFillArray(x, ret)
    ret
  }
  
  def jacksonKernel(order: Int): Array[R] = {
    import math._
    val Mp = order+1d
    Array.tabulate(order) { m => 
      (1/Mp)*((Mp-m)*cos(Pi*m/Mp) + sin(Pi*m/Mp)/tan(Pi/Mp))
    }
  }
  
  // Coefficients c_m of linear combination of moments for weight calculation
  //   F = \sum mu_m c_m = \int rho(e) fn(e)
  def expansionCoefficients(order: Int, de: R, fn: R => R): Array[R] = {
    import math._
    val kernel = jacksonKernel(order)
    val ret = Array.fill[R](order)(0)
    for (e <- -1.0+de to 1.0-de by de) {
      val t = KPM.chebyshevArray(e, order)
      val f = fn(e) / (Pi * sqrt(1 - e*e))
      for (m <- 0 until order) {
        ret(m) += (if (m == 0) 1 else 2) * kernel(m) * t(m) * f * de
      }
    }
    // for ((f, i) <- ret.zipWithIndex) println("Coeff %d = %g".format(i, f))
    ret
  }

  // --- Plotting stuff ---
    
  // converts moments into a density of states
  def densityOfStates(moments: Array[R], kernel: Array[R], e: R): R = {
    import math._
    val order = moments.size
    val t = KPM.chebyshevArray(e, order)
    var ret = 0d
    for (m <- 0 until order) {
      ret += (if (m == 0) 1 else 2) * kernel(m) * moments(m) * t(m)
    }
    ret / (Pi * sqrt(1 - e*e))
  }

  def eigenvaluesApprox(order: Int, range: Array[R], kpm: KPM): Array[R] = {
    // val mu = time("Calculating %d moments (slow)".format(order))(kpm.momentsExact())
    val r = kpm.randomVector()
    val kernel = jacksonKernel(order)
    val mu = kpm.momentsStochastic(order, r)
    range.map(densityOfStates(mu, kernel, _))
  }
  
  def eigenvaluesExact(H: PackedSparse[S]): Array[Double] = {
    val (v, w) = H.toDense.map(_.toComplexd).eig
    v.map(_.re).toArray.sorted
  }
  
    // \int dx x^{moment} y(x) 
  def integrate(xs: Array[R], ys: Array[R], moment: Int): Array[R] = {
    require(xs.size >= 2 && xs.size == ys.size)
    require((1 until xs.size).forall(i => xs(i-1) < xs(i)))
    xs.indices.toArray.foldMapLeft(0: R) { (acc, i) =>
      val x0 = if (i > 0)         xs(i-1) else xs(i)
      val x2 = if (i < xs.size-1) xs(i+1) else xs(i)
      val x = 0.5*(x0 + x2)
      val dx = 0.5 * (x2 - x0)
      acc + dx*math.pow(x, moment)*ys(i)
    }
  }
  
  // \int dx x^{moment} \sum \delta(x - x_i)
  def integrateDeltas(xs: Array[R], deltas: Array[Double], moment: Int): Array[R] = {
    require(xs.size >= 2)
    require((1 until xs.size).forall(i => xs(i-1) < xs(i)))
    require((1 until deltas.size).forall(i => deltas(i-1) <= deltas(i)))
    require(deltas.head >= xs.head && deltas.last <= xs.last)
    val ret = new Array[R](xs.length)
    var j = 0
    var acc = 0d
    for (i <- xs.indices) {
      val binEnd: R = if (i < xs.size-1) avg(xs(i), xs(i+1)) else xs(i)
      while (j < deltas.size && deltas(j) <= binEnd) {
        acc += math.pow(deltas(j), moment)
        j += 1
      }
      ret(i) = acc
    }
    ret
  }
  

  def range(npts: Int): Array[R] = {
    Array.tabulate[R](npts) { i =>
      2.0 * (i+0.5) / npts - 1.0
    }
  }
  
  def mkPlot(name: String): Plot = {
    val plot = new Plot(name)
    scikit.util.Utilities.frame(plot.getComponent(), plot.getTitle())
    plot
  }
  
  def plotHistogram(plot: Plot, a: Array[R]) {
    val hist = new scikit.dataset.Histogram(0.01)
    a.foreach { hist.accum(_) }
    plot.registerBars("Density", hist, java.awt.Color.BLACK)
  }
  
  def plotLines(plot: Plot, data: (Array[R], Array[R]), name: String = "data", color: java.awt.Color = java.awt.Color.BLACK) {
    val pts = new scikit.dataset.PointSet(data._1.map(_.toDouble), data._2.map(_.toDouble))
    plot.registerLines(name, pts, color)
  }

}


abstract class GenKPM(val H: PackedSparse[S], val nrand: Int, val seed: Int) {
  val rand = new util.Random(seed)
  val n = H.numRows

  def momentsStochastic(order: Int, r: Dense[S]): Array[R]
  def functionAndGradient(r: Dense[S], c: Array[R], grad: PackedSparse[S]): R
  
  def momentsExact(order: Int): Array[R] = {
    val ret = Array.fill[R](order)(0)
    val r = dense(n, nrand)
    for (i <- 0 until n by nrand) {
      r.transform(_ => 0)
      for (j <- 0 until nrand; if i+j < n)
        r(i+j, j) = 1
      ret += momentsStochastic(order, r)
    }
    ret.transform(_ * nrand)     // undo division in momentsStochastic
    ret(0) = n                   // first two moments are handled specially
    ret(1) = H.trace.re
    ret
  }
  
  def momentsExactDense(order: Int): Array[R] = {
    val ret = Array.fill[R](order)(0)
    
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
  
  def gradientExactDense(c: Array[R], grad: PackedSparse[S]) {
    val order = c.size
    grad.clear()
    val Ht = H.tran
    
    val u0 = eye(n).toDense  // U_{0}
    val u1 = Ht.toDense * 2  // U_{1}
    val u2 = dense(n, n)     // U_{2}

    for ((i, j) <- grad.definedIndices) {
      grad(i, j) += 1 * c(1) * u0(i, j)
    }

    for (m <- 2 until order) {
      // u0 = U_{m-2}(H)
      // u1 = U_{m-1}(H)
      // u2 = U_m(H) = 2 H t1 - t0
      u2 := u0
      u2.gemm(2, Ht, u1, -1)
      
      for ((i, j) <- grad.definedIndices) {
        grad(i, j) += m * c(m) * u1(i, j)
      }

      u0 := u1
      u1 := u2
    }
  }

  // Vector with randomly values, uniform among (1, i, -1, -i)
  def randomVector(): Dense[S] = {
    val r  = dense(n, nrand)
    val a = Array[S#A](1, I, -1, -I)
    r.fill(a(rand.nextInt(4)))
  }
  
  def randomGaussianVector(): Dense[S] = {
    val r  = dense(n, nrand)
    r.fill((rand.nextGaussian() + I*rand.nextGaussian) / math.sqrt(2))
    // normalizing random vectors may slightly improve accuracy
//    for (j <- 0 until r.numCols) {
//      val c = r(::, j).norm2.re
//      for (i <- 0 until r.numRows) {
//        r(i, j) *= math.sqrt(n / c) 
//      }
//    }
//    println(r(::,0).norm2.re)
    r
  }
}


// Kernel polynomial method

class KPM(H: PackedSparse[S], nrand: Int, seed: Int = 0) extends GenKPM(H, nrand, seed) {
  
  override def momentsStochastic(order: Int, r: Dense[S]): Array[R] = momentsStochasticAux(order, r)._1
  
  // Returns: (mu(m), alpha_{M-2}, alpha_{M-1})
  def momentsStochasticAux(order: Int, r: Dense[S]): (Array[R], Dense[S], Dense[S]) = {
    val mu = Array.fill[R](order)(0)
    mu(0) = n                   // Tr[T_0[H]] = Tr[1]
    mu(1) = H.trace.re          // Tr[T_1[H]] = Tr[H]
    
    val a0 = dense(n, nrand)
    val a1 = dense(n, nrand)
    val a2 = dense(n, nrand)
    
    a0 := r                      // T_0[H] |r> = 1 |r>
    a1 :=* (H, r)                // T_1[H] |r> = H |r>
    
    for (m <- 2 to order-1) {
      a2 := a0; a2.gemm(2, H, a1, -1)  // alpha_m = T_m[H] r = 2 H a1 - a0

      mu(m) = (r dagDot a2).re / nrand
      a0 := a1
      a1 := a2
    }
    
    (mu, a0, a1)
  }
  
  override def functionAndGradient(r: Dense[S], c: Array[R], grad: PackedSparse[S]): R = {
    val order = c.size
    grad.clear()
    
    val a2 = dense(n, nrand)
    val (mu, a0, a1) = momentsStochasticAux(order, r)
    
    val b2 = dense(n, nrand)
    val b1 = dense(n, nrand)
    val b0 = r * c(order - 1)
    
    // need special logic since (mu_1) is calculated exactly
    for (i <- 0 until grad.numRows) { grad(i, i) += c(1) }
    def cp(m: Int): R = if (m == 1) 0 else c(m)
    
    // cache defined indices for speed
    val (indicesI, indicesJ) = {
      val (i, j) = grad.definedIndices.unzip
      (i.toArray, j.toArray)
    }
    
    for (m <- order-2 to 0 by -1) {
      // a0 = alpha_{m}
      // b0 = beta_{m}

      if (nrand > 1) {
        for ((i, j) <- grad.definedIndices; k <- 0 until nrand) {
          grad(i, j) += (if (m == 0) 1 else 2) * b0(i, k).conj * a0(j, k) / nrand
        }
      }
      // equivalent to above, but much faster. b2 is used as a temporary vector.
      else {
        if (m == 0) (b2 := b0) else (b2 :=* (2, b0))
        for (iter <- 0 until indicesI.length) {
          grad.scalar.maddTo(true, b2.data, indicesI(iter), a0.data, indicesJ(iter), grad.data, iter)
        }
      }
      
      a2 := a1
      b2 := b1
      a1 := a0
      b1 := b0
      a0 := a2; a0.gemm(2, H, a1, -1)                   // a0 = 2 H a1 - a2 
      b0 :=* (cp(m), r); b0.gemm(2, H, b1, 1); b0 -= b2 // b0 = c(m) r + 2 H b1 - b2
    }
    
    (c, mu).zipped.map(_*_).sum
  }
}
