package kip.projects.quantum.kpm

import smatrix._
import scala.util.Random

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

object KPMUtil {
  def jacksonKernel(order: Int): Array[Double] = {
    import math._
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
  
  def expansionCoefficients(M: Int, quadPts: Int, f: Double => Double, es: EnergyScale): Array[Double] = {
    // TODO: replace with DCT-II, f -> fp
    var fp = Array.fill[Double](quadPts)(0d)
    val T = Array.fill[Double](M)(0d)
    for (i <- 0 until quadPts) {
      val x_i = math.cos(math.Pi * (i+0.5) / quadPts)
      chebyshevFillArray(x_i, T)
      for (m <- 0 until M) {
        fp(m) += T(m) * f(es.unscale(x_i))
      }
    }
    
    val kernel = jacksonKernel(M)
    var ret = Array.fill[Double](M)(0)
    for (m <- 0 until M) {
      ret(m) = /* es.mag * */ (if (m == 0) 1 else 2) * kernel(m) * fp(m) / quadPts 
    }
    ret
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


trait ComplexKPM {
  type Cd = Scalar.ComplexDbl
  def moments(M: Int, r: Dense[Cd], H: PackedSparse[Cd], es: EnergyScale): Array[Double]
  def functionAndGradient(c: Array[Double], r: Dense[Cd], H: PackedSparse[Cd], es: EnergyScale): (Double, PackedSparse[Cd])
}

object ComplexKPMCpu extends ComplexKPM {
  import Constructors.complexDbl._
  
  // Returns: (mu(m), alpha_{M-2}, alpha_{M-1})
  def momentsAux(M: Int, r: Dense[Cd], Hs: PackedSparse[Cd]): (Array[Double], Dense[Cd], Dense[Cd]) = {
    val n = r.numRows
    val s = r.numCols
    
    val mu = Array.fill[Double](M)(0)
    mu(0) = n                   // Tr[T_0[H]] = Tr[1]
    mu(1) = Hs.trace.re          // Tr[T_1[H]] = Tr[H]
    
    val a0 = dense(n, s)
    val a1 = dense(n, s)
    val a2 = dense(n, s)
    
    a0 := r                      // T_0[H] |r> = 1 |r>
    a1 :=* (Hs, r)                // T_1[H] |r> = H |r>
    
    for (m <- 2 to M-1) {
      a2 := a0; a2.gemm(2, Hs, a1, -1)  // alpha_m = T_m[H] r = 2 H a1 - a0
      
      mu(m) = (r dagDot a2).re
      a0 := a1
      a1 := a2
    }
    
    (mu, a0, a1)
  }
  
  def moments(M: Int, r: Dense[Cd], H: PackedSparse[Cd], es: EnergyScale): Array[Double] = {
    momentsAux(M, r, es.scale(H))._1
  }
  
  def functionAndGradient(c: Array[Double], r: Dense[Cd], H: PackedSparse[Cd], es: EnergyScale): (Double, PackedSparse[Cd]) = {
    val n = r.numRows
    val s = r.numCols
    val M = c.size
    
    val Hs = es.scale(H)
    val dE_dHs = Hs.duplicate.clear()
    
    val a2 = dense(n, s)
    val (mu, a0, a1) = momentsAux(M, r, Hs)
    
    val b2 = dense(n, s)
    val b1 = dense(n, s)
    val b0 = r * c(M - 1)
    
    // need special logic since `mu_1` is calculated exactly
    // note that `mu_0` does not contribute to derivative
    for (i <- 0 until n) { dE_dHs(i, i) += c(1) }
    def cp(m: Int): Double = if (m == 1) 0 else c(m)
    
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
    
    val E = (c, mu).zipped.map(_*_).sum
    val dE_dH = dE_dHs.transform(_ / es.mag)
    (E, dE_dH)
  }
}
