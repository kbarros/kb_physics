package kip.projects.quantum.qmd

import Units._
import kip.math.Vec3
import kip.projects.quantum.kpm.ComplexKPM
import kip.projects.quantum.kpm.KPMUtil
import smatrix._
import smatrix.Constructors.complexDbl._
import kip.projects.quantum.kpm.EnergyScale


class TbHamiltonian(pot: Potential, lat: Lattice, x: Array[Vec3]) {
  val nOrbs = pot.numOrbitalsPerSite
  val nAtoms = lat.numAtoms
  val n = nOrbs * nAtoms
  val nspin = 2.0

  val H = {
    val H     = sparse(n, n)
    val h     = Array.ofDim[Double](nOrbs, nOrbs)
    val tmp   = Array.ofDim[Double](nOrbs, nOrbs)
    for (i <- 0 until lat.numAtoms;
         j <- lat.neighbors(i, x, pot.rcut);
         if (j <= i)) {
      pot.fillTBHoppings(lat.displacement(x(i), x(j)), h, tmp, tmp, tmp)
      for (o1 <- 0 until nOrbs;
           o2 <- 0 until nOrbs) {
        H(i*nOrbs+o1, j*nOrbs+o2) = h(o1)(o2)
        H(j*nOrbs+o2, i*nOrbs+o1) = h(o1)(o2)
      }
    }
    
    // val hermitDev = math.sqrt((q.matrix - q.matrix.dag).norm2.abs)
    // require(hermitDev < 1e-6, "Found non-hermitian hamiltonian! Deviation: "+hermitDev)
    
    H.toPacked
  }
  
  // returns dX_dr for every atomic position r
  def chainGradient(dX_dH: PackedSparse[Scalar.ComplexDbl]): Array[Vec3] = {
    val ret = Array.fill(nAtoms)(Vec3.zero)
    val tmp   = Array.ofDim[Double](nOrbs, nOrbs)
    val dh_dx = Array.ofDim[Double](nOrbs, nOrbs)
    val dh_dy = Array.ofDim[Double](nOrbs, nOrbs)
    val dh_dz = Array.ofDim[Double](nOrbs, nOrbs)
    for (i <- 0 until lat.numAtoms;
         j <- lat.neighbors(i, x, pot.rcut);
         if (j <= i)) {
      pot.fillTBHoppings(lat.displacement(x(i), x(j)), tmp, dh_dx, dh_dy, dh_dz)
      for (o1 <- 0 until nOrbs;
           o2 <- 0 until nOrbs) {
        val dX_dH_ij = (dX_dH(i*nOrbs+o1, j*nOrbs+o2) + dX_dH(j*nOrbs+o2, i*nOrbs+o1)).re
        // force applied by atom i on atom j
        val dX_dr_j = Vec3(dh_dx(o1)(o2), dh_dy(o1)(o2), dh_dz(o1)(o2)) * dX_dH_ij
        ret(j) += dX_dr_j
        ret(i) -= dX_dr_j
      }
    }
    ret
  }
  
  def localFermiEnergy(x: Double, T: Double, mu: Double) = {
    import math.{abs, exp, log}
    val alpha = (x-mu)/(kB*abs(T))
    if (T == 0.0 || abs(alpha) > 20) {
      if (x < mu) nspin*(x-mu) else 0.0
    }
    else {
      -kB*T*nspin*log(1 + exp(-alpha))
    }
  }
  
  def localFermiDensity(x: Double, T: Double, mu: Double) = {
    import math.{abs, exp, log}
    val alpha = (x-mu)/(kB*abs(T))
    if (T == 0.0 || abs(alpha) > 20) {
      if (x < mu) nspin else 0.0
    }
    else {
      // println(s"at x=$x got n=${1.0/(exp(alpha)+1.0)}")
      1.0/(exp(alpha)+1.0)
    }
  }
  
  def energyAndForce(kpm: ComplexKPM, fd: ComplexKPM.ForwardData, mu: Double, T: Double): (Double, Array[Vec3]) = {
    var E = 0.0
    val f = Array.fill(nAtoms)(Vec3.zero)
    
    // electronic part
    E += kpm.function(fd, localFermiEnergy(_, T, mu))
    val dE_dr = chainGradient(kpm.gradient(fd, localFermiEnergy(_, T, mu)))
    for (i <- 0 until lat.numAtoms)
      f(i) -= dE_dr(i)
    
    // pair part
    for (i <- 0 until lat.numAtoms;
         j <- lat.neighbors(i, x, pot.rcut);
         if (j < i)) {
      val del = lat.displacement(x(i), x(j))
      val r = del.norm
      E += pot.phi(r)
      val f_ij = del * (-pot.dphi_dr(r) / r)
      f(j) += f_ij
      f(i) -= f_ij
    }
    
    (E, f)
  }
  
  def fillingFraction(kpm: ComplexKPM, fd: ComplexKPM.ForwardData, mu: Double, T: Double) = {
    kpm.function(fd, localFermiDensity(_, T, mu)) / (nspin*n)
  }
  
  def fillingGradient(kpm: ComplexKPM, fd: ComplexKPM.ForwardData, mu: Double, T: Double) = {
    chainGradient(kpm.gradient(fd, localFermiDensity(_, T, mu)))
  }
}
