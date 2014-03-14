package kip.projects.quantum.qmd2

import kip.math.Vec3
import kip.projects.quantum.qmd.{Potential, GoodwinSi, Lattice}
import kip.projects.quantum.qmd.Units._
import kip.projects.quantum.kpm2.KPMComplex
import kip.projects.quantum.kpm2.KPMUtil
import kip.projects.quantum.kpm2.EnergyScale
import kip.projects.quantum.kpm2.SparseCsrComplex
import kip.projects.quantum.kpm2.SparseCooComplex


class TbHamiltonian(pot: Potential, lat: Lattice, x: Array[Vec3]) {
  val nOrbs = pot.numOrbitalsPerSite
  val nAtoms = lat.numAtoms
  val n = nOrbs * nAtoms
  val nspin = 2.0
  
  val H = new SparseCsrComplex(n, n)
  val H_builder = new SparseCooComplex(n, n)
  
  def buildHamiltonian() {
    H_builder.clear()
    val h     = Array.ofDim[Double](nOrbs, nOrbs)
    val tmp   = Array.ofDim[Double](nOrbs, nOrbs)
    for (i <- 0 until lat.numAtoms;
         j <- lat.neighbors(i, x, pot.rcut);
         if (j <= i)) {
      pot.fillTBHoppings(lat.displacement(x(i), x(j)), h, tmp, tmp, tmp)
      for (o1 <- 0 until nOrbs;
           o2 <- 0 until nOrbs) {
        H_builder.add(i*nOrbs+o1, j*nOrbs+o2, h(o1)(o2), 0.0)
        if (i != j)
          H_builder.add(j*nOrbs+o2, i*nOrbs+o1, h(o1)(o2), 0.0)
      }
    }
    H.fromCoo(H_builder)
  }
  
  // returns dX_dr for every atomic position r
  def chainGradient(dX_dH: SparseCsrComplex): Array[Vec3] = {
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
        val dX_dH_ij = dX_dH.get_re(i*nOrbs+o1, j*nOrbs+o2)
        val dX_dH_ji = if (i == j) 0.0 else dX_dH.get_re(j*nOrbs+o2, i*nOrbs+o1)
        val dX_dr_j = Vec3(dh_dx(o1)(o2), dh_dy(o1)(o2), dh_dz(o1)(o2)) * (dX_dH_ij + dX_dH_ji)
        ret(j) += dX_dr_j
        ret(i) -= dX_dr_j
      }
    }
    ret
  }
  
  def grouping(i: Int, s: Int): Int = {
    val norb = pot.numOrbitalsPerSite
    norb*lat.grouping(i/norb, x, s/norb) + (i%norb)
  }
  
  def localFermiEnergy(x: Double, T: Double, mu: Double) = {
    import math.{abs, exp, log}
    val alpha = (x-mu)/(kB*abs(T))
    if (true || T == 0.0 || abs(alpha) > 20) {
      if (x < mu) nspin*(x-mu) else 0.0
    }
    else {
      -kB*T*nspin*log(1 + exp(-alpha))
    }
  }
  
  def localFermiDensity(x: Double, T: Double, mu: Double) = {
    import math.{abs, exp, log}
    val alpha = (x-mu)/(kB*abs(T))
    if (true || T == 0.0 || abs(alpha) > 20) {
      if (x < mu) nspin else 0.0
    }
    else {
      1.0/(exp(alpha)+1.0)
    }
  }
  
  def energy(kpm: KPMComplex, mu: Double, T: Double): Double = {
    // electronic part
    var E = kpm.eval(localFermiEnergy(_, T, mu))
    for (i <- 0 until lat.numAtoms;
         j <- lat.neighbors(i, x, pot.rcut);
         if (j < i)) {
      val del = lat.displacement(x(i), x(j))
      val r = del.norm
      E += pot.phi(r)
    }
    E
  }
  
  // mu and fillingFraction must correspond
  def energyAtFixedFilling(kpm: KPMComplex, mu: Double, fillingFraction: Double, T: Double): Double = {
    val E = energy(kpm, mu, T)
    val occupiedStates = fillingFraction*nspin*n
    E + mu*occupiedStates
  }
  
  def force(kpm: KPMComplex, mu: Double, T: Double): Array[Vec3] = {
    // electronic part
    val f = Array.fill(nAtoms)(Vec3.zero)
    val dE_dr = chainGradient(kpm.gradient(localFermiEnergy(_, T, mu)))
    for (i <- 0 until lat.numAtoms)
      f(i) -= dE_dr(i)
    // pair part
    for (i <- 0 until lat.numAtoms;
         j <- lat.neighbors(i, x, pot.rcut);
         if (j < i)) {
      val del = lat.displacement(x(i), x(j))
      val r = del.norm
      val f_ij = del * (-pot.dphi_dr(r) / r)
      f(j) += f_ij
      f(i) -= f_ij
    }
    f
  }
  
  // TODO: generalize to Fermi function
  // TODO: search for band gap
  def findChemicalPotential(kpm: KPMComplex, fillingFraction: Double): Double = {
    val (x, irho) = KPMUtil.integratedDensityFunction(kpm.gamma, kpm.es)
    val i = irho.indexWhere(_ > fillingFraction*n)
    require(0.0 < fillingFraction && fillingFraction < 1.0, s"Filling fraction $fillingFraction out of range (-1, 1)")
    require(0 < i && i < irho.size-1, s"Could not find mu for filling fraction $fillingFraction")
    // println(s"mu1=${x(i-1)} n1=${irho(i-1)/8} mu2=${x(i)} n2=${irho(i)/8}")
    (x(i) + x(i-1)) / 2.0
  }
  
  def fillingFraction(kpm: KPMComplex, mu: Double, T: Double) = {
    kpm.eval(localFermiDensity(_, T, mu)) / (nspin*n)
  }
}
