package kip.projects.quantum.qmd

import smatrix._
import Constructors.complexDbl._
import kip.math.Vec3
import kip.javasim.Random

object TbHamiltonian extends App {
  testDimer()
  
  def testDimer() {
    val pot = new GoodwinSi(rcut=10)
    val r = pot.r0
    
    val lat = new LinearChain(numAtoms=2, spacing=r)
    
    val delta = Vec3(1.1, 2.3, -1.8).normalize * (-r)
    
    def buildHamiltonian(delta: Vec3) = {
      new TbHamiltonian(pot, lat, Array(Vec3.zero, delta))
    }
    val tbh = buildHamiltonian(delta)
    
    def sqr(x: Double) = x*x
    val eig1 = Array(
      pot.Δsp + pot.hppπ,
      pot.Δsp + pot.hppπ,
      pot.Δsp - pot.hppπ,
      pot.Δsp - pot.hppπ,
      (pot.hssσ - pot.hppσ)/2 + pot.Δsp/2 + math.sqrt(sqr(pot.Δsp/2-(pot.hssσ+pot.hppσ)/2)+sqr(pot.hspσ)),
      (pot.hssσ - pot.hppσ)/2 + pot.Δsp/2 - math.sqrt(sqr(pot.Δsp/2-(pot.hssσ+pot.hppσ)/2)+sqr(pot.hspσ)),
      -(pot.hssσ - pot.hppσ)/2 + pot.Δsp/2 + math.sqrt(sqr(pot.Δsp/2+(pot.hssσ+pot.hppσ)/2)+sqr(pot.hspσ)),
      -(pot.hssσ - pot.hppσ)/2 + pot.Δsp/2 - math.sqrt(sqr(pot.Δsp/2+(pot.hssσ+pot.hppσ)/2)+sqr(pot.hspσ))
    ).sorted
    val eig2 = tbh.H.toDense.eig._1.toArray.map(_.re).sorted
    
    val eigDev = math.sqrt((eig1 zip eig2).map(e => sqr(e._1-e._2)).sum)
    println(s"Eigenvalue deviation: $eigDev")
    
    val rand = new Random(0)
    val rmat = tbh.H.duplicate.transform(x => if (x == 0) 0 else rand.nextGaussian())
    val de_dh = (rmat + rmat.dag).toPacked
    
    // fake energy with simple, random dependence on elements of H 
    def energy(H: PackedSparse[Scalar.ComplexDbl]) = (de_dh.toDense dagDot H.toDense).re
    
    val dx = 0.0001
    def deriv(dir: Vec3): Double = {
      val ep = energy(buildHamiltonian(delta + dir*dx).H)
      val em = energy(buildHamiltonian(delta - dir*dx).H)
      (ep - em) / (2*dx)
    }
    val forceDiscrete = -Vec3(deriv(Vec3(1,0,0)), deriv(Vec3(0,1,0)), deriv(Vec3(0,0,1)))
    val forceAnalytical = tbh.forces(de_dh)(1)
    
    // println(s"$forceDiscrete $forceAnalytical")
    println(s"Force deviation: ${forceAnalytical-forceDiscrete}")
  }
}

class TbHamiltonian(pot: Potential, lat: Lattice, pos: Array[Vec3]) {
  val nOrbs = pot.numOrbitalsPerSite
  val nAtoms = lat.numAtoms
  val n = nOrbs * nAtoms
  
  val H = {
    val H     = sparse(n, n)
    val h     = Array.ofDim[Double](nOrbs, nOrbs)
    val tmp   = Array.ofDim[Double](nOrbs, nOrbs)
    for (i <- 0 until lat.numAtoms;
         j <- lat.neighbors(i, pos, pot.rcut);
         if (j <= i)) {
      pot.fillTBHoppings(lat.displacement(pos(i), pos(j)), h, tmp, tmp, tmp)
      for (o1 <- 0 until nOrbs;
           o2 <- 0 until nOrbs) {
        H(i*nOrbs+o1, j*nOrbs+o2) = h(o1)(o2)
        H(j*nOrbs+o2, i*nOrbs+o1) = h(o1)(o2)
      }
    }
    H.toPacked
  }
  
  
  // returns forces (-dE/dr_i) on every atom
  
  def forces(dE_dH: PackedSparse[Scalar.ComplexDbl]): Array[Vec3] = {
    val f = Array.fill(nAtoms)(Vec3.zero)
    
    val tmp   = Array.ofDim[Double](nOrbs, nOrbs)
    val dh_dx = Array.ofDim[Double](nOrbs, nOrbs)
    val dh_dy = Array.ofDim[Double](nOrbs, nOrbs)
    val dh_dz = Array.ofDim[Double](nOrbs, nOrbs)
    
    for (i <- 0 until lat.numAtoms;
         j <- lat.neighbors(i, pos, pot.rcut);
         if (j <= i)) {
      pot.fillTBHoppings(lat.displacement(pos(i), pos(j)), tmp, dh_dx, dh_dy, dh_dz)
      for (o1 <- 0 until nOrbs;
           o2 <- 0 until nOrbs) {
        val dE_dH_ij = (dE_dH(i*nOrbs+o1, j*nOrbs+o2) + dE_dH(j*nOrbs+o2, i*nOrbs+o1)).re
        // force applied by atom i on atom j
        val f_ij = - Vec3(dh_dx(o1)(o2), dh_dy(o1)(o2), dh_dz(o1)(o2)) * dE_dH_ij
        f(j) += f_ij
        f(i) -= f_ij
      }
    }
    f
  }
}
