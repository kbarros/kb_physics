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
      println("Num interactions = "+is.size)
      (for (i <- is; j <- i.compatibleInteractions(is)) yield i.cutoff(j)).max
  }

  def potentialEnergy(): Double = {
    var ret = 0.0
    val cutoff = globalCutoff()
    volume.buildCells(atomsPerCell, atoms)
    
    for (a1 <- atoms) {
      ret += a1.potential1(this)

      for (a2 <- volume.atomsInRange(a1, cutoff)) {
	if (a1 != a2) {
	  ret += a1.potential2(this, a2)
        }
      }
    }
    ret
  }
  
  def kineticEnergy(): Double = {
    var ret = 0.0
    for (a <- atoms) {
      ret += 0.5*a.mass*(sqr(a.vx) + sqr(a.vy) + sqr(a.vz))
    }
    ret
  }
  
  def calculateForces() {
    for (a <- atoms) {
      a.fx = 0
      a.fy = 0
      a.fz = 0
    }
    val cutoff = globalCutoff()
    volume.buildCells(atomsPerCell, atoms)

    for (a1 <- atoms) {
      
      val f = a1.force1(this)
      a1.fx += f.x
      a1.fy += f.y
      a1.fz += f.z
      
      for (a2 <- volume.atomsInRange(a1, cutoff)) {
        if (a1 != a2) {
	  val (f1, f2) = a1.force2(this, a2)
          a1.fx += f1.x
          a1.fy += f1.y
          a1.fz += f1.z
          a2.fx += f2.x
          a2.fy += f2.y
          a2.fz += f2.z
        }
      }
    }
  }
  
  
  def step() {
    println("atoms" + atoms(0) + " " + atoms(1))
    println("pot: " + potentialEnergy())
    integrator.singleStep(this)
    time += integrator.dt
  }
  
  def visualize(viz: Visualizer, radius: Atom => Double, color: Atom => java.awt.Color) {
    viz.setBounds(volume.bounds)
    viz.setParticles(atoms.map(a => Visualizer.Sphere(a.pos, radius(a), color(a))))
  }
}
