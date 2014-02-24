package kip.projects.quantum.qmd

import kip.math.Vec3
import smatrix._
import smatrix.Constructors.complexDbl._


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
