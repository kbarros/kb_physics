package kip.md

import scala.math._
import kip.math._
import kip.math.Math._
import kip.graphics._


class World(var volume: Volume, var atoms: Seq[Atom], var integrator: Integrator) {
  var atomsPerCell = 4
  var time = 0.0
  volume.wrapAtoms(atoms)
  
  def globalCutoff() = {
      val is = atoms.toSet[Atom].flatMap(_.tag.inter2)
      (for (i <- is; j <- i.compatibleInteractions(is)) yield i.cutoff(j)).max
  }

  def potentialEnergy(): Double = {
    var ret = 0.0
    val cutoff = globalCutoff()
    volume.buildCells(atomsPerCell, atoms)
    
    for (a1 <- atoms) {
      for (i <- a1.tag.inter1) {
        ret += i.potential(this, a1)
      }
      
      for (a2 <- volume.atomsInRange(a1, cutoff)) {
	if (a1.idx < a2.idx) {
          for (i1 <- a1.tag.inter2;
               i2 <- i1.compatibleInteractions(a2.tag.inter2)) {
            ret += i1.potential(this, a1, i2, a2)
          }
        }
      }
    }
    ret
  }
  
  def kineticEnergy(): Double = {
    var ret = 0.0
    for (a <- atoms) {
      ret += 0.5*a.mass*(sqr(a.v.x) + sqr(a.v.y) + sqr(a.v.z))
    }
    ret
  }
  
  def temperature(): Double = {
    val dof = 3.0 * atoms.size
    (2.0/dof) * kineticEnergy()
  }
  
  def forceOnObject(inter1: Interaction1): Vec3 = {
    val f = Vec3.zero.toMutable
    for (a <- atoms)
      inter1.accumForce(this, a, f)
    f
  }
  
  def calculateForces() {
    for (a <- atoms)
      a.f.reset()
    val cutoff = globalCutoff()
    volume.buildCells(atomsPerCell, atoms)
    val f = mutable.Vec3.zero
    
    for (a1 <- atoms) {
      for (i1 <- a1.tag.inter1) {
        i1.accumForce(this, a1, a1.f)
      }
      
      for (a2 <- volume.atomsInRange(a1, cutoff)) {
        if (a1 != a2) {
          for (i1 <- a1.tag.inter2;
               i2 <- i1.compatibleInteractions(a2.tag.inter2)) {
            if (!i1.symmetric) {
              i1.accumForce(this, a1, i2, a2, a1.f)
            }
            else if (a1.idx < a2.idx) {
              f.reset()
              i1.accumForce(this, a1, i2, a2, f)
              a1.f += f
              a2.f -= f
            }
          }
        }
      }
    }
  }
  
  
  def step(steps: Int) {
    for (i <- 0 until steps) {
      integrator.singleStep(this)
      time += integrator.dt
    }
  }
}
