package kip.projects.quantum.kpm

import scala.util.Random
import math._
import smatrix._

class EnergyScale(val lo: Double, val hi: Double) {
  val avg = (hi + lo) / 2.0
  val mag = (hi - lo) / 2.0 
  
  def scale(x: Double) = {
    (x - avg) / mag 
  }
  
  def scale(x: PackedSparse[Scalar.ComplexDbl]) = {
    val ret = x.duplicate
    for (i <- 0 until ret.numRows) {
      ret(i, i) -= avg
    }
    ret.transform(_ / mag)
    ret
  }
  
  def unscale(x: Double) = {
    x * mag + avg
  }
}

// TODO: merge KPMUtil and KPM?

object KPMUtil {
  def jacksonKernel(order: Int): Array[Double] = {
    val Mp = order+1d
    Array.tabulate(order) { m => 
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
  def expansionCoefficients(M: Int, quadPts: Int, f: Double => Double, es: EnergyScale): Array[Double] = {
    // TODO: replace with DCT-II, f -> fp
    var fp = Array.fill[Double](quadPts)(0d)
    val T = Array.fill[Double](M)(0d)
    for (i <- 0 until quadPts) {
      val x_i = math.cos(math.Pi * (i+0.5) / quadPts)
      chebyshevFillArray(x_i, T)
      for (m <- 0 until M) {
        fp(m) += f(es.unscale(x_i)) * T(m) 
      }
    }
    
    val kernel = jacksonKernel(M)
    var ret = Array.fill[Double](M)(0)
    for (m <- 0 until M) {
      ret(m) = /* es.mag * */ (if (m == 0) 1 else 2) * kernel(m) * fp(m) / quadPts 
    }
    ret
  }
  
  // Transformation of moments mu suitable for reconstructing density of states
  def momentTransform(mu: Array[Double], quadPts: Int): Array[Double] = {
    val M = mu.size
    val T = new Array[Double](M)
    val mup = (mu, jacksonKernel(M)).zipped.map(_*_)
    val gamma = new Array[Double](quadPts)
    
    // TODO: replace with DCT-III, mup -> gamma
    for (i <- 0 until quadPts) {
      val x_i = cos(Pi * (i+0.5) / quadPts)
      chebyshevFillArray(x_i, T) // T_m(x_i) = cos(m pi (i+1/2) / quadPts)
      for (m <- 0 until M) {
        gamma(i) += (if (m == 0) 1 else 2) * mup(m) * T(m)
      }
    }
    gamma
  }
  
  // Estimate \int_lo^hi rho(x) g(x) dx 
  def densityProduct(gamma: Array[Double], g: Double => Double, es: EnergyScale): Double = {
    val quadPts = gamma.size
    var ret = 0.0
    for (i <- 0 until quadPts) {
      val x_i = cos(Pi * (i+0.5) / quadPts)
      ret += gamma(i) * g(es.unscale(x_i))
    }
    ret / quadPts
  }
  
  // Returns density of states rho(x) at Chebyshev points x
  def densityFunction(gamma: Array[Double], es: EnergyScale): (Array[Double], Array[Double]) = {
    val quadPts = gamma.size
    val x = new Array[Double](quadPts)
    val rho = new Array[Double](quadPts)
    for (i <- quadPts-1 to 0 by -1) {
      val x_i = cos(Pi * (i+0.5) / quadPts)
      x(quadPts-1-i) = es.unscale(x_i)
      rho(quadPts-1-i) = gamma(i) / (Pi * sqrt(1-x_i*x_i) * es.mag)
    }
    (x, rho)
  }
  
  // Returns density of states \int theta(x-x') rho(x') dx' at Chebyshev points x
  def integratedDensityFunction(gamma: Array[Double], es: EnergyScale): (Array[Double], Array[Double]) = {
    val quadPts = gamma.size
    val x = new Array[Double](quadPts)
    val irho = new Array[Double](quadPts)
    var acc = 0.0
    for (i <- quadPts-1 to 0 by -1) {
      val x_i = cos(Pi * (i+0.5) / quadPts)
      x(quadPts-1-i) = es.unscale(x_i)
      irho(quadPts-1-i) = (acc+0.5*gamma(i)) / quadPts
      acc += gamma(i) 
    }
    (x, irho)
  }

  // Vectors with random elements in {1, i, -1, -i}
  def uncorrelatedVectors(n: Int, s: Int, rand: Random): Dense[Scalar.ComplexDbl] = {
    import Constructors.complexDbl._
    val r = dense(n, s)
    val a = Array[Complexd](1, I, -1, -I).map(_ / math.sqrt(s))
    r.fill(a(rand.nextInt(4)))
  }
  
  // Nearly orthogonal vectors with mostly zeros
  def correlatedVectors(n: Int, s: Int, grouping: Int => Int, rand: Random): Dense[Scalar.ComplexDbl] = {
    import Constructors.complexDbl._
    val r = dense(n, s)
    val a = Array[Complexd](1, I, -1, -I)
    r.tabulate { (i, j) =>
      val g = grouping(i)
      require(0 <= g && g < s)
      if (g == j) a(rand.nextInt(4)) else 0
    }
  }
  
  def allVectors(n: Int): Dense[Scalar.ComplexDbl] = {
    import Constructors.complexDbl._
    val r = dense(n, n)
    for (i <- 0 until n) r(i, i) = 1
    r
  }
  
  def energyScale(H: PackedSparse[Scalar.ComplexDbl]) = {
    val eig_min = H.eig(nev=1, which="SR", tol=1e-4)._1.apply(0).re
    val eig_max = H.eig(nev=1, which="LR", tol=1e-4)._1.apply(0).re
    val slack = 0.01 * (eig_max - eig_min)
    new EnergyScale(lo = eig_min - slack, hi = eig_max + slack)
  }
}

object ComplexKPM {
  type Cd = Scalar.ComplexDbl
  val quadAccuracy = 10
  case class ForwardData(Hs: PackedSparse[Cd], es: EnergyScale, r: Dense[Cd], mu: Array[Double], gamma: Array[Double], aM2: Dense[Cd], aM1: Dense[Cd])
}

trait ComplexKPM {
  import ComplexKPM._
  
  def forward(M: Int, r: Dense[Cd], H: PackedSparse[Cd], es: EnergyScale): ForwardData
  def reverse(fd: ForwardData, coeff: Array[Double]): PackedSparse[Cd]
  
  def function(fd: ForwardData, f: Double=>Double): Double = {
    KPMUtil.densityProduct(fd.gamma, f, fd.es)
  }

  def gradient(fd: ForwardData, f: Double=>Double): PackedSparse[Cd] = {
    val M = fd.mu.size
    val quadPts = M * quadAccuracy
    val coeff = KPMUtil.expansionCoefficients(M, quadPts, f, fd.es)
    reverse(fd, coeff)
  }
}

object ComplexKPMCpu extends ComplexKPM {
  import ComplexKPM._
  import Constructors.complexDbl._
  
  def forward(M: Int, r: Dense[Cd], H: PackedSparse[Cd], es: EnergyScale): ForwardData = {
    val n = r.numRows
    val s = r.numCols
    val Hs = es.scale(H)
    
    val mu = Array.fill[Double](M)(0)
    mu(0) = n                     // Tr[T_0[H]] = Tr[1]
    mu(1) = Hs.trace.re           // Tr[T_1[H]] = Tr[H]
    
    val a0 = dense(n, s)
    val a1 = dense(n, s)
    val a2 = dense(n, s)
    a0 := r                       // T_0[H] |r> = 1 |r>
    a1 :=* (Hs, r)                // T_1[H] |r> = H |r>
    
    for (m <- 2 to M-1) {
      a2 := a0; a2.gemm(2, Hs, a1, -1)  // alpha_m = T_m[H] r = 2 H a1 - a0
      
      mu(m) = (r dagDot a2).re
      a0 := a1
      a1 := a2
    }
    
    val gamma = KPMUtil.momentTransform(mu, quadAccuracy*M)
    ForwardData(Hs, es, r, mu, gamma, aM2=a0, aM1=a1)
  }
  
  def reverse(fd: ForwardData, coeff: Array[Double]): PackedSparse[Cd] = {
    import fd._
    val n = r.numRows
    val s = r.numCols
    val M = coeff.size

    val dE_dHs = Hs.duplicate.clear()
    
    val a2 = dense(n, s)
    val a1 = aM1.duplicate
    val a0 = aM2.duplicate
    
    val b2 = dense(n, s)
    val b1 = dense(n, s)
    val b0 = r * coeff(M - 1)
    
    // need special logic since `mu_1` is calculated exactly
    // note that `mu_0` does not contribute to derivative
    for (i <- 0 until n) { dE_dHs(i, i) += coeff(1) }
    def cp(m: Int): Double = if (m == 1) 0 else coeff(m)
    
    // cache defined indices for speed
    val (indicesI, indicesJ) = {
      val (i, j) = dE_dHs.definedIndices.unzip
      (i.toArray, j.toArray)
    }
    
    for (m <- M-2 to 0 by -1) {
      // a0 = alpha_{m}
      // b0 = beta_{m}
      
      if (s > 1) {
        for ((i, j) <- dE_dHs.definedIndices; k <- 0 until s) {
          dE_dHs(i, j) += (if (m == 0) 1 else 2) * b0(i, k).conj * a0(j, k)
        }
      }
      // equivalent to above, but much faster. b2 is used as a temporary vector.
      else {
        if (m == 0) (b2 := b0) else (b2 :=* (2, b0))
        for (iter <- 0 until indicesI.length) {
          dE_dHs.scalar.maddTo(true, b2.data, indicesI(iter), a0.data, indicesJ(iter), dE_dHs.data, iter)
        }
      }
      
      a2 := a1
      b2 := b1
      a1 := a0
      b1 := b0
      a0 := a2; a0.gemm(2, Hs, a1, -1)                   // a0 = 2 H a1 - a2 
      b0 :=* (cp(m), r); b0.gemm(2, Hs, b1, 1); b0 -= b2 // b0 = c(m) r + 2 H b1 - b2
    }
    
    dE_dHs.transform(_ / es.mag)
  }
}
