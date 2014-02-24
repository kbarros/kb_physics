package kip.projects.quantum.qmd

import Units._
import math._
import scala.util.Random
import kip.math.Vec3
import kip.projects.quantum.kpm.ComplexKPM
import kip.projects.quantum.kpm.KPMUtil


class TbMD(x: Array[Vec3], v: Array[Vec3], pot: Potential, lat: Lattice, kpm: ComplexKPM, M: Int) {
  val N = lat.numAtoms

  def filling(mu: Double) = {    
    def fillingFn(x: Double) = {
      if (x < mu) 1.0 else 0.0
    }
    val tbh = new TbHamiltonian(pot, lat, x)
    val es = KPMUtil.energyScale(tbh.H)
    val c = KPMUtil.expansionCoefficients(M=M, quadPts=4*M, f=fillingFn, es=es)
    val r = KPMUtil.allVectors(tbh.n)
    kpm.functionAndGradient(c, r, tbh.H, es)._1 / tbh.n 
  }
  
  def energyAndForces(mu: Double, T: Double): (Double, Array[Vec3]) = {
    def energyFn(x: Double) = {
      val nspin = 2.0
      val alpha = (x-mu)/(kB*max(T,+0.0))
      if      (alpha < -20) nspin*(x - mu)
      else if (alpha > +20) 0.0
      else                 -kB*T*nspin*log(1 + exp(-alpha))
    }
    val tbh = new TbHamiltonian(pot, lat, x)
    val es = KPMUtil.energyScale(tbh.H)
    val c = KPMUtil.expansionCoefficients(M=M, quadPts=4*M, f=energyFn, es=es)
    val r = KPMUtil.allVectors(tbh.n)
    val (e, de_dH) = kpm.functionAndGradient(c, r, tbh.H, es)
    // val f = tbh.forces(de_dH) // strange Scala bug
    (e, tbh.forces(de_dH))
  }
  
  def timestep(f: Array[Vec3], m: Double, gamma: Double, T: Double, dt: Double, rand: Random) {
    import math.sqrt
    
    // update momentum at fixed position
    for (i <- 0 until N) {
      val eta = Vec3(rand.nextGaussian(), rand.nextGaussian(), rand.nextGaussian())
      v(i) += f(i)*(-dt/m) + v(i)*(-gamma*dt) + eta * sqrt(2*kB*T*gamma*dt/m)
    }
    
    // update position at fixed momentum
    for (i <- 0 until N) {
      x(i) += v(i) * dt
    }
  }
}

