package kip.projects.quantum.kpm2

import scala.util.Random
import math._

class EnergyScale(val lo: Double, val hi: Double) {
  val avg = (hi + lo) / 2.0
  val mag = (hi - lo) / 2.0
  
  override def toString() = s"EnergyScale(lo=$lo, hi=$hi)"
  
  def scale(x: Double) = {
    (x - avg) / mag 
  }
  
  def scale(x: SparseCsrComplex): Unit = {
    x += (-avg, 0)
    x *= (1.0/mag, 0)
  }
  
  def unscale(x: Double) = {
    x * mag + avg
  }
}


object KPMUtil {
  def jacksonKernel(M: Int): Array[Double] = {
    val Mp = M+1.0
    Array.tabulate(M) { m => 
      (1/Mp)*((Mp-m)*cos(Pi*m/Mp) + sin(Pi*m/Mp)/tan(Pi/Mp))
    }
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
  
  // Returns array c such that:
  //    \int_lo^hi rho(x) f(x) = \sum mu_m c_m
  def expansionCoefficients(M: Int, Mq: Int, f: Double => Double, es: EnergyScale): Array[Double] = {
    // TODO: replace with DCT-II, f -> fp
    var fp = Array.fill[Double](Mq)(0d)
    val T = Array.fill[Double](M)(0d)
    for (i <- 0 until Mq) {
      val x_i = math.cos(math.Pi * (i+0.5) / Mq)
      chebyshevFillArray(x_i, T)
      for (m <- 0 until M) {
        fp(m) += f(es.unscale(x_i)) * T(m) 
      }
    }
    
    val kernel = jacksonKernel(M)
    var ret = Array.fill[Double](M)(0)
    for (m <- 0 until M) {
      ret(m) = (if (m == 0) 1 else 2) * kernel(m) * fp(m) / Mq 
    }
    ret
  }
  
  // Transformation of moments mu suitable for reconstructing density of states
  def momentTransform(mu: Array[Double], Mq: Int): Array[Double] = {
    val M = mu.size
    val T = new Array[Double](M)
    val mup = (mu, jacksonKernel(M)).zipped.map(_*_)
    val gamma = new Array[Double](Mq)
    
    // TODO: replace with DCT-III, mup -> gamma
    for (i <- 0 until Mq) {
      val x_i = cos(Pi * (i+0.5) / Mq)
      chebyshevFillArray(x_i, T) // T_m(x_i) = cos(m pi (i+1/2) / Mq)
      for (m <- 0 until M) {
        gamma(i) += (if (m == 0) 1 else 2) * mup(m) * T(m)
      }
    }
    gamma
  }
  
  // Estimate \int_lo^hi rho(x) g(x) dx 
  def densityProduct(gamma: Array[Double], g: Double => Double, es: EnergyScale): Double = {
    val Mq = gamma.size
    var ret = 0.0
    for (i <- 0 until Mq) {
      val x_i = cos(Pi * (i+0.5) / Mq)
      ret += gamma(i) * g(es.unscale(x_i))
    }
    ret / Mq
  }
  
  // Returns density of states rho(x) at Chebyshev points x
  def densityFunction(gamma: Array[Double], es: EnergyScale): (Array[Double], Array[Double]) = {
    val Mq = gamma.size
    val x = new Array[Double](Mq)
    val rho = new Array[Double](Mq)
    for (i <- Mq-1 to 0 by -1) {
      val x_i = cos(Pi * (i+0.5) / Mq)
      x(Mq-1-i) = es.unscale(x_i)
      rho(Mq-1-i) = gamma(i) / (Pi * sqrt(1-x_i*x_i) * es.mag)
    }
    (x, rho)
  }
  
  // Returns density of states \int theta(x-x') rho(x') dx' at Chebyshev points x
  def integratedDensityFunction(gamma: Array[Double], es: EnergyScale): (Array[Double], Array[Double]) = {
    val Mq = gamma.size
    val x = new Array[Double](Mq)
    val irho = new Array[Double](Mq)
    var acc = 0.0
    for (i <- Mq-1 to 0 by -1) {
      val x_i = cos(Pi * (i+0.5) / Mq)
      x(Mq-1-i) = es.unscale(x_i)
      irho(Mq-1-i) = (acc+0.5*gamma(i)) / Mq
      acc += gamma(i) 
    }
    (x, irho)
  }
  
  def energyScale(H: SparseCsrComplex) = {
    val Hp = H.toSmatrix()
    val eig_min = Hp.eig(nev=1, which="SR", tol=1e-4)._1.apply(0).re
    val eig_max = Hp.eig(nev=1, which="LR", tol=1e-4)._1.apply(0).re
    val slack = 0.01 * (eig_max - eig_min)
    new EnergyScale(lo = eig_min - slack, hi = eig_max + slack)
  }
}


abstract class KPMComplex(val H: SparseCsrComplex, val s: Int, val M: Int, val Mq: Int) {
  val n = H.numRows
  val Hs = new SparseCsrComplex(n, n)
  val dX_dH = new SparseCsrComplex(n, n)
  val R = new DenseComplex(n, s)
  val mu = new Array[Double](M)
  var gamma = new Array[Double](Mq)
  var es: EnergyScale = null
  
  def forward(es: EnergyScale)
  def gradient(f: Double=>Double): SparseCsrComplex  
  
  // Vectors with random elements in {1, i, -1, -i}
  def uncorrelatedVectors(rand: Random) {
    val x = 1.0 / math.sqrt(s)
    for (j <- 0 until R.numCols;
         i <- 0 until R.numRows) {
      rand.nextInt(4) match {
        case 0 => R.set(i, j, +x,  0)
        case 1 => R.set(i, j,  0, +x)
        case 2 => R.set(i, j, -x,  0)
        case 3 => R.set(i, j,  0, -x)
      }
    }
  }
  
  // Nearly orthogonal vectors with mostly zeros
  def correlatedVectors(grouping: Int => Int, rand: Random) {
    R.zero()
    for (i <- 0 until R.numRows) {
      val g = grouping(i)
      require(0 <= g && g < R.numCols)
      rand.nextInt(4) match {
        case 0 => R.set(i, g, +1,  0)
        case 1 => R.set(i, g, -1,  0)
        case 2 => R.set(i, g,  0, +1)
        case 3 => R.set(i, g,  0, -1)
      }
    }
  }
  
  def allVectors(): Unit = {
    require(R.numRows == R.numCols)
    R.zero()
    for (i <- 0 until R.numRows) {
      R.set(i, i, 1, 0)
    }
  }
  
  def eval(f: Double=>Double): Double = {
    KPMUtil.densityProduct(gamma, f, es)
  }
}


class KPMComplexCpu(H: SparseCsrComplex, s: Int, M: Int, Mq: Int) extends KPMComplex(H, s, M, Mq) {
  val alphaM2 = new DenseComplex(n, s)
  val alphaM1 = new DenseComplex(n, s)

  def forward(es: EnergyScale) {
    this.es = es
    Hs.fromCsr(H)
    es.scale(Hs)
    
    import smatrix._
    import Constructors.complexDbl._
    
    val _R = R.toSmatrix()
    val _Hs = Hs.toSmatrix()
    
    mu(0) = n                      // Tr[T_0[H]] = Tr[1]
    mu(1) = _Hs.trace.re           // Tr[T_1[H]] = Tr[H]
    
    val a0 = dense(n, s)
    val a1 = dense(n, s)
    val a2 = dense(n, s)
    a0 := _R                       // T_0[H] |r> = 1 |r>
    a1 :=* (_Hs, _R)               // T_1[H] |r> = H |r>
    
    for (m <- 2 to M-1) {
      a2 := a0; a2.gemm(2, _Hs, a1, -1)  // alpha_m = T_m[H] r = 2 H a1 - a0
      
      mu(m) = (_R dagDot a2).re
      a0 := a1
      a1 := a2
    }
    
    gamma = KPMUtil.momentTransform(mu, Mq)
    alphaM2.fromSmatrix(a0)
    alphaM1.fromSmatrix(a1)
  }
  
  def gradient(f: Double=>Double): SparseCsrComplex = {
    val coeff = KPMUtil.expansionCoefficients(M, Mq, f, es)
    
    import smatrix._
    import Constructors.complexDbl._
    
    val _R = R.toSmatrix
    val _Hs = Hs.toSmatrix
    val _dX_dHs = _Hs.duplicate.clear()
    
    val a2 = dense(n, s)
    val a1 = alphaM1.toSmatrix
    val a0 = alphaM2.toSmatrix
    
    val b2 = dense(n, s)
    val b1 = dense(n, s)
    val b0 = _R * coeff(M - 1)
    
    // need special logic since `mu_1` is calculated exactly
    // note that `mu_0` does not contribute to derivative
    for (i <- 0 until n) { _dX_dHs(i, i) += coeff(1) }
    def cp(m: Int): Double = if (m == 1) 0 else coeff(m)
    
    for (m <- M-2 to 0 by -1) {
      // a0 = alpha_{m}
      // b0 = beta_{m}
      for ((i, j) <- _dX_dHs.definedIndices; k <- 0 until s) {
        _dX_dHs(i, j) += (if (m == 0) 1 else 2) * b0(i, k).conj * a0(j, k)
      }      
      a2 := a1
      b2 := b1
      a1 := a0
      b1 := b0
      a0 := a2; a0.gemm(2, _Hs, a1, -1)                    // a0 = 2 H a1 - a2 
      b0 :=* (cp(m), _R); b0.gemm(2, _Hs, b1, 1); b0 -= b2 // b0 = c(m) r + 2 H b1 - b2
    }
    
    dX_dH.fromSmatrix(_dX_dHs.transform(_ / es.mag))
    dX_dH
  }
}
