package kip.projects.bridging

import math._
import scala.util.Random

case class Vector(x: Double, y: Double) {
  val abs2 = x*x + y*y
  val abs = sqrt(abs2)
  def *(a: Double) = Vector(x*a, y*a)
  def +(v: Vector) = Vector(x+v.x, y+v.y)
  def -(v: Vector) = Vector(x-v.x, y-v.y)
  def dot(v: Vector) = x*v.x + y*v.y
}


class SubMD2d(val ncols: Int, val nrows: Int, val a: Double, val dt: Double, val m: Double = 1, val sigma: Double = 1, val rc: Double = 3) {
  val rnd = new Random()
  val n = ncols*nrows // number of atoms
  
  // general Strain type would be a full 3x3 matrix. for now, only handle stretch along x direction
  type Strain = Double
  type Stress = Double
  var strain: Strain = 1
  
  // hexagonal crystal in equilibrium
  //
  //   o   o   o   o
  // o   o   o   o 
  //   o   o   o   o
  // o   o   o   o 
  // <-a->
  
  // width/height of reference box
  val w0 = a*ncols
  val h0 = (a*sqrt(3)/2)*nrows

  // physical (x, y) coordinates in deformed box (i.e., box is transformed according to strain tensor) 
  val p = Array.fill(n)(Vector(0, 0))
  initializeCrystalPositions()
  
  // physical velocity in deformed box (i.e., box is transformed according to strain tensor) 
  val v = Array.fill(n)(Vector(0, 0))
  
  // apply new strain to all particles (position and velocities)
  def applyStrain(newStrain: Strain) {
    val oldStrain = this.strain
    for (i <- 0 until n) {
      p(i) = forwardTransform(newStrain, backwardTransform(oldStrain, p(i)))
      v(i) = forwardTransform(newStrain, backwardTransform(oldStrain, v(i)))
    }
    this.strain = newStrain
  }
  
  // volume of deformed box 
  def volume() = {
    val v0 = w0*h0
    val jacobian = strain // more generally: jacobian of deformation gradient
    v0 * jacobian
  }
  
  // transform vector from reference to physical coordinates (according to strain)
  def forwardTransform(strain: Strain, r: Vector) = Vector(r.x*strain, r.y)
  
  // transform vector from physical to reference coordinates (according to strain)
  def backwardTransform(strain: Strain, r: Vector) = Vector(r.x/strain, r.y)
  
  // initialize particles to strained hexagonal crystal
  def initializeCrystalPositions() {
    require (nrows % 2 == 0)
    for (iy <- 0 until nrows; ix <- 0 until ncols) {
      val i = iy*ncols + ix
      // crystal position in reference coordinates 
      val p0 = Vector(x = ix*a + (iy%2)*(a/2), y = iy*(a*sqrt(3)/2))
      // transform to physical coordinates
      p(i) = forwardTransform(strain, p0)
    }
  }
  
  def sqr(x: Double) = x*x
  def cube(x: Double) = x*x*x
  def mod(x: Double, y: Double) = ((x % y) + y) % y      // mathematical mod; result in range [0, y)
  def wrap(x: Double, y: Double) = mod(x + y/2, y) - y/2 // wrap x into range [-y/2, y/2)

  def lennardJonesPotential(r: Double): Double = {
    val sr = sigma / r
    val sr2 = sqr(sr)
    val sr6 = cube(sr2)
    4 * (sqr(sr6) - sr6)
  }
  
  val lennardJonesShift = lennardJonesPotential(rc)
  
  // Shifted truncated Lennard Jones potential:
  // if r < rc
  //   4 ((s/r)^12 - (s/r)^6) - 4 ((s/rc)^12 - (s/rc)^6)
  // else
  //   0
  def lennardJonesPotentialTruncated(r: Double) = {
    if (r < rc) lennardJonesPotential(r) - lennardJonesShift
    else 0
  }
  
  // Truncated Lennard Jones force:
  //
  // if r < rc
  //   (\vec r) (24/s^2) (2 (s/r)^14 - (s/r)^8)
  // else
  //   0
  def lennardJonesForceTruncated(vr: Vector): Vector = {
    val r = vr.abs
    if (r < rc) {
      val sr = sigma / r
      val sr2 = sqr(sr)
      val sr6 = cube(sr2)
      val sr12 = sqr(sr6)
      vr * ( (24/sqr(sigma)) * (2 * sr12 - sr6) * sr2 )
    }
    else Vector(0, 0)
  }
  
  def kineticEnergy() = (m/2) * v.map(_.abs2).sum
  
  
  def wrapAtomPosition(i: Int) {
    // particle position in reference (undeformed) coordinates
    val p0 = backwardTransform(strain, p(i))
    // wrap atom position in reference cubic volume, then re-apply deformation
    p(i) = forwardTransform(strain, Vector(mod(p0.x, w0), mod(p0.y, h0)))
  } 

  // Calculates displacement vector (p0(i) - p0(j)) in reference coordinates
  def referenceDisplacement(i: Int, j: Int): Vector = {
    // particle positions in reference crystal
    val p0i = backwardTransform(strain, p(i))
    val p0j = backwardTransform(strain, p(j))
    
    // use minimal image in cubic volume with periodic boundary conditions
    Vector(x = wrap(p0i.x - p0j.x, w0), y = wrap(p0i.y - p0j.y, h0))    
  }
  
  // Calculates displacement vector (p(i) - p(j)) in physical coordinates.
  // Assumes that the nearest image in the reference volume is also the
  // nearest image in the deformed volume.
  def displacement(i: Int, j: Int): Vector = {
    forwardTransform(strain, referenceDisplacement(i, j))
  }
  
  // naive summation over all pairs
  def potentialEnergy() = {
    var ret = 0d
    for (i <- 0 until n; j <- (i+1) until n) {
      val r = displacement(i, j).abs
      ret += lennardJonesPotentialTruncated(r)
    }
    ret
  }
  
  def energy() = kineticEnergy() + potentialEnergy()
  
  def netForce(i: Int) = {
    var ret = Vector(0, 0)
    for (j <- 0 until n; if j != i) {
      ret += lennardJonesForceTruncated(displacement(i, j))
    }
    ret
  }
  
  def verletStep() {
    // update velocities based on forces
    for (i <- 0 until n) {
      val f = netForce(i)
      v(i) += f * (dt / m)
    }
    
    // update positions based on velocity
    for (i <- 0 until n) {
      p(i) += v(i) * dt
      wrapAtomPosition(i)
    }
  }
  
  // re-initialize particle velocities to achieve target kinetic energy
  def applyKineticEnergy(targetKinetic: Double) {
    // randomize velocities
    v.transform(_ => Vector(rnd.nextGaussian(), rnd.nextGaussian()))
    
    // adjust velocities so that net momentum is zero
    var pnet = Vector(0, 0)
    for (i <- 0 until n) pnet += v(i) * m
    for (i <- 0 until n) v(i) -= pnet * (1.0 / (m * n))
    
    // rescale velocities to match target kinetic energy
    val ke = kineticEnergy()
    val scale = sqrt(targetKinetic / kineticEnergy())
    v.transform(_*scale)
  }
  
  // the virial stress is consistent with the macroscale Cauchy stress
  //   http://en.wikipedia.org/wiki/Virial_stress
  // for now, compute only (1,1) component
  def virialStress(): Stress = {
    var ret: Stress = 0d
    // kinetic part
    for (i <- 0 until n) {
      ret += - m * v(i).x * v(i).x
    }
    // potential part
    for (i <- 0 until n; j <- i+1 until n) {
      val r = displacement(i, j)
      val f = lennardJonesForceTruncated(r) // force on atom i due to atom j
      ret += - r.x * f.x // just need the (1,1) component of stress tensor 
    }
    ret / volume()
  }
  
  // see definitions here: http://en.wikipedia.org/wiki/Stress_(mechanics)
  // 
  // the 1st Piola-Kirchoff stress tensor is: P = (det F) F^-1 sigma
  // where sigma is the Cauchy stress and F is the deformation gradient dx / dx^0
  //
  // note that sigma is symmetric, but P is not, because it connects reference
  // to material coordinates (ie., it is a "two-point" tensor).
  //
  // in the present case of our 1d problem, where F = diag(a, 1, 1), this is a no-op.
  // 
  def convertCauchyToFirstPiolaKirchoff(s: Stress): Stress = {
    s
  }
}