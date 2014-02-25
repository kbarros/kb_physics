package kip.projects.quantum.qmd

import Units._
import math._
import scala.util.Random
import kip.math.Vec3
import kip.projects.quantum.kpm.ComplexKPM
import kip.projects.quantum.kpm.KPMUtil


class TbMD(x: Array[Vec3], v: Array[Vec3], pot: Potential, lat: Lattice, kpm: ComplexKPM, M: Int) {
  
  def timestep(f: Array[Vec3], m: Double, gamma: Double, T: Double, dt: Double, rand: Random) {
    import math.sqrt
    
    // update momentum at fixed position
    for (i <- 0 until lat.numAtoms) {
      val eta = Vec3(rand.nextGaussian(), rand.nextGaussian(), rand.nextGaussian())
      v(i) += f(i)*(-dt/m) + v(i)*(-gamma*dt) + eta * sqrt(2*kB*T*gamma*dt/m)
    }
    
    // update position at fixed momentum
    for (i <- 0 until lat.numAtoms) {
      x(i) += v(i) * dt
    }
  }
}

