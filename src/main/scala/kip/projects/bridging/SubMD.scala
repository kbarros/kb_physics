package kip.projects.bridging

import math._

class SubMD2d(val ncols: Int, val nrows: Int, val a: Double, val dt: Double, val m: Double = 1, val sigma: Double = 1, val rc: Double = 3) {
  val n = ncols*nrows // number of atoms
  
  case class Vector(x: Double, y: Double) {
    val abs2 = x*x + y*y 
    val abs = sqrt(abs2)
    def *(a: Double) = Vector(x*a, y*a)
    def +(v: Vector) = Vector(x+v.x, y+v.y)
  }
  
  // physical (x, y) coordinates in deformed box (i.e., box is transformed according to strain tensor) 
  val p = new Array[Vector](n)

  // physical velocity in deformed box (i.e., box is transformed according to strain tensor) 
  val v = new Array[Vector](n)
  
  // hexagonal crystal in equilibrium
  //
  //   o   o   o   o
  // o   o   o   o 
  //   o   o   o   o
  // o   o   o   o 
  // <-a->
  
  // width/height of undeformed box
  val wu = a*ncols + 0.321
  val hu = (a*sqrt(3)/2)*nrows
  
  // initialize particles in hexagonal crystal
  require (nrows % 2 == 0)
  for (iy <- 0 until nrows; ix <- 0 until ncols) {
    val i = iy*ncols + ix
    p(i) = Vector(x = ix*a + (iy%2)*(a/2), y = iy*(a*sqrt(3)/2))
    v(i) = Vector(0, 0)
  }
  
  // general Strain type would be a full 3x3 matrix. for now, only handle stretch along x direction
  type Strain = Double
  var strain: Strain = 1
  
  // apply new strain to all particles (position and velocities)
  def applyStrain(newStrain: Strain) {
    val oldStrain = this.strain
    for (i <- 0 until n) {
      p(i) = forwardTransform(newStrain, backwardTransform(oldStrain, p(i)))
      v(i) = forwardTransform(newStrain, backwardTransform(oldStrain, v(i)))
    }
    this.strain = newStrain
  }
  
  // apply strain tensor to vector
  def forwardTransform(strain: Strain, r: Vector) = Vector(r.x*strain, r.y)
  
  // apply inverse strain tensor to vector
  def backwardTransform(strain: Strain, r: Vector) = Vector(r.x/strain, r.y)
  
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
    // backtransform position to cubic volume of reference crystal 
    val pu = backwardTransform(strain, p(i))
    // wrap atom position in cubic volume, and then reapply strain deformation
    p(i) = forwardTransform(strain, Vector(mod(pu.x, wu), mod(pu.y, hu)))
  } 

  // Calculates displacement (p(i) - p(j)) accounting for periodic boundary conditions
  // of deformed simulation volume. Uses the heuristic that the nearest image in the
  // reference volume is also the nearest image in the deformed volume.
  def atomDisplacement(i: Int, j: Int): Vector = {
    // backtransform positions to cubic volume of reference crystal 
    val pi = backwardTransform(strain, p(i))
    val pj = backwardTransform(strain, p(j))
    
    // find distance in reference crystal, according to periodic boundary conditions
    val r = Vector(x = wrap(pi.x - pj.x, wu), y = wrap(pi.y - pj.y, hu))
    
    // reapply strain deformation
    forwardTransform(strain, r)
  }
  
  // naive summation over all pairs
  def potentialEnergy() = {
    var ret = 0d
    for (i <- 0 until n; j <- (i+1) until n) {
      val r = atomDisplacement(i, j).abs
      ret += lennardJonesPotentialTruncated(r)
    }
    ret
  }
  
  def energy() = kineticEnergy() + potentialEnergy()
  
  def netForce(i: Int) = {
    var ret = Vector(0, 0)
    for (j <- 0 until n; if j != i) {
      ret += lennardJonesForceTruncated(atomDisplacement(i, j))
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
}