package kip.projects.quantum

import kip.util.Util.{time, notime}
import kip.math.Math._
import kip.enrich._
import smatrix._
import ctor._


object KPM {
  import scikit.graphics.dim2._
  
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
    val pts = new scikit.dataset.PointSet(data._1, data._2)
    plot.registerLines(name, pts, color)
  }
  
    
  def eigenvaluesExact(H: PackedSparse[S]): Array[R] = {
    val (v, w) = H.toDense.eig
    v.map(_.re).toArray.sorted
  }

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
  
  // \int dx x^{moment} \sum \delta(x - x_i)
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

class KPM(val H: PackedSparse[S], val order: Int, val nrand: Int, val seed: Int = 0) {
  val rand = new util.Random(seed)
  val n = H.numRows
  
  lazy val jacksonKernel: Array[R] = {
    import math._
    val Mp = order+1d
    Array.tabulate(order) { m => 
      (1/Mp)*((Mp-m)*cos(Pi*m/Mp) + sin(Pi*m/Mp)/tan(Pi/Mp))
    }
  }
  lazy val kernel = jacksonKernel
  
  // Coefficients c_m of linear combination of moments for weight calculation
  //   F = \sum mu_m c_m = \int rho(e) fn(e)
  def expansionCoefficients(de: Double, fn: Double => Double): Array[Double] = {
    import math._
    val ret = Array.fill[Double](order)(0d)
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
  
  // converts moments into a density of states
  def densityOfStates(mu: Array[R], e: R): R = {
    import math._
    val t = KPM.chebyshevArray(e, order)
    var ret = 0d
    for (m <- 0 until order) {
      ret += (if (m == 0) 1 else 2) * kernel(m) * mu(m) * t(m)
    }
    ret / (Pi * sqrt(1 - e*e))
  }
  
  def momentsExact(): Array[R] = {
    val ret = Array.fill(order)(0d)
    
    // TODO: test with sparse
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
  
  // Vector with randomly values, uniform among (1, i, -1, -i)
  def randomVector(): Dense[S] = {
    val r  = dense(n, nrand)
    r.fill(Seq[S#A](1, I, -1, -I).apply(rand.nextInt(4)))
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
  
  // Returns: (mu(m), alpha_{M-2}, alpha_{M-1})
  def momentsStochastic(r: Dense[S]): (Array[R], Dense[S], Dense[S]) = {
    val mu = Array.fill(order)(0d)
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
  
  def gradientExact(c: Array[R], grad: PackedSparse[S]) {
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
  
  def functionAndGradient(r: Dense[S], c: Array[R], grad: PackedSparse[S]): R = {
    grad.clear()
    
    val a2 = dense(n, nrand)
    val (mu, a0, a1) = momentsStochastic(r)
    
    val b2 = dense(n, nrand)
    val b1 = dense(n, nrand)
    val b0 = r * c(order - 1)
    
    // need special logic since (mu_1) is calculated exactly
    for (i <- 0 until grad.numRows) { grad(i, i) += c(1) }
    def cp(m: Int) = if (m == 1) 0d else c(m)
    
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
      // equivalent to above, but much faster. b3 is used as a temporary vector.
      else {
        println("fail; need to cplx conjugate b1 and iterate over nrand")
//        if (m > 1) (b3 :=* (2, b1)) else (b3 := b1)
//        for (iter <- 0 until indicesI.length) {
//          grad.scalar.maddTo(false, b3.data, indicesI(iter), a0.data, indicesJ(iter), grad.data, iter)
//        }
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
  
  val range: Array[R] = {
    val nx = 5*order 
    Array.tabulate[R](nx) { i =>
      2.0 * (i+0.5) / nx - 1.0
    }
  }
  
  def eigenvaluesApprox(kernel: Array[R]): Array[R] = {
//    // Test reordering of sums is correct
//    val r1 = randomVector()
//    val forward = momentsStochastic(r1)
//    val rho = range.map(densityOfStates(forward._1, _))
//    println("version A " + ((rho, range).zipped.map(_*_).sum * (range(1) - range(0))))
//    println("version B " + functionAndGradient(e => e, null)) // different random vector => different result

//    val mu = time("Calculating %d moments (slow)".format(order))(momentsExact())
    val r = randomVector()
    val (mu, _, _) = notime("Calculating %d moments (stoch) of N=%d matrix".format(order, H.numRows))(momentsStochastic(r))
    range.map(densityOfStates(mu, _))
  }
}
