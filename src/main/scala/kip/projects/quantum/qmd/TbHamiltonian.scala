package kip.projects.quantum.qmd

import Units._
import kip.math.Vec3
import kip.projects.quantum.kpm.ComplexKPM
import kip.projects.quantum.kpm.KPMUtil
import smatrix._
import smatrix.Constructors.complexDbl._


class TbHamiltonian(pot: Potential, lat: Lattice, x: Array[Vec3]) {
  val nOrbs = pot.numOrbitalsPerSite
  val nAtoms = lat.numAtoms
  val n = nOrbs * nAtoms
  
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
    H.toPacked
  }
  
  
  def energyAndForces(mu: Double, T: Double, kpm: ComplexKPM, M: Int): (Double, Array[Vec3]) = {
    import math.{max, exp, log}
    
    var e = 0.0
    val f = Array.fill(nAtoms)(Vec3.zero)

    // electronic part
    def energyFn(x: Double) = {
      val nspin = 2.0
      val alpha = (x-mu)/(kB*max(T,+0.0))
      if (alpha < -20)
        nspin*(x - mu)
      else if (alpha > +20)
        0.0
      else
        -kB*T*nspin*log(1 + exp(-alpha))
    }
    val es = KPMUtil.energyScale(H)
    val c = KPMUtil.expansionCoefficients(M, 4*M, energyFn, es)
    val r = KPMUtil.allVectors(n)
    val (electronicEnergy, de_dH) = kpm.functionAndGradient(c, r, H, es)
    e += electronicEnergy
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
        val de_dH_ij = (de_dH(i*nOrbs+o1, j*nOrbs+o2) + de_dH(j*nOrbs+o2, i*nOrbs+o1)).re
        // force applied by atom i on atom j
        val f_ij = - Vec3(dh_dx(o1)(o2), dh_dy(o1)(o2), dh_dz(o1)(o2)) * de_dH_ij
        f(j) += f_ij
        f(i) -= f_ij
      }
    }
    
    // pair part
    for (i <- 0 until lat.numAtoms;
         j <- lat.neighbors(i, x, pot.rcut);
         if (j < i)) {
      val del = lat.displacement(x(i), x(j))
      val r = del.norm
      e += pot.phi(r)
      val f_ij = del * (-pot.dphi_dr(r) / r)
      f(j) += f_ij
      f(i) -= f_ij
    }

    (e, f)
  }
  
  def filling(mu: Double, kpm: ComplexKPM, M: Int) = {
    def fillingFn(x: Double) = {
      if (x < mu) 1.0 else 0.0
    }
    val es = KPMUtil.energyScale(H)
    val c = KPMUtil.expansionCoefficients(M, 4*M, fillingFn, es)
    val r = KPMUtil.allVectors(n)
    kpm.functionAndGradient(c, r, H, es)._1 / n
  }

}
