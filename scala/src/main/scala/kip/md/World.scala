package kip.md

import scala.math._
import kip.math.Math._


class World {
  var periodic = true
  var temperature = 1.0
  var integrator = new Verlet(this, 0.1)

  private var _atoms = Seq[Atom]()
  private var _grid: PointGrid2d[Atom] = _
  private var _globalCutoff: Double = _
  
  def setSize(Lx: Double, Ly: Double) {
    val atomsPerCell = 4
    val cols = sqrt(_atoms.size / atomsPerCell).toInt
    _grid = new PointGrid2d(Lx, cols, periodic)
  }
  
  def atoms = _atoms
  def atoms_= (atoms: Seq[Atom]) {
    _atoms = atoms
    _globalCutoff = {
      val is = atoms.toSet[Atom].flatMap(_.tag.inter2)
      println("Num interactions = "+is.size)
      (for (i <- is; j <- i.compatibleInteractions(is)) yield i.cutoff(j)).max
    }
  }
  
  def wrapAtoms() {
    def wrap(x: Double, wx: Int): (Double, Int) = {
      val L = _grid.L
      if (x >= 2*L || x < -L) {
        println("Simulation exploded")
        exit(-1)
      }
      if (x >= _grid.L) (x-L, wx+1)
      else if (x < 0) (x+L, wx-1)
      else (x, wx)
    }
    
    for (a <- _atoms) {
      val (x, wx) = wrap(a.x, a.wx)
      val (y, wy) = wrap(a.y, a.wy)
      val (z, wz) = wrap(a.z, a.wz)
      a.x = x
      a.y = y
      a.z = z
      a.wx = wx
      a.wy = wy
      a.wz = wz
    }
  }
  
  def potentialEnergy: Double = {
    var ret = 0.0
    _grid.loadPoints(_atoms)
    for (a1 <- _atoms) {
      ret += a1.potential1

      for (a2 <- _grid.pointOffsetsWithinRange(a1, _globalCutoff)) {
	if (a1 != a2)
	  ret += a1.potential2(a2)
      }
    }
    ret
  }
  
  def kineticEnergy: Double = {
    var ret = 0.0
    for (a <- _atoms) {
      ret += 0.5*a.mass*(sqr(a.vx) + sqr(a.vy) + sqr(a.vz))
    }
    ret
  }
  
  def calculateForces() {
    for (a <- _atoms) {
      a.fx = 0
      a.fy = 0
      a.fz = 0
    }
    
    _grid.loadPoints(_atoms)
    for (a1 <- _atoms) {
      
      val f = a1.force1
      a1.fx += f.x
      a1.fy += f.y
      a1.fz += f.z
      
      for (a2 <- _grid.pointOffsetsWithinRange(a1, _globalCutoff)) {
        if (a1 != a2) {
	  val (f1, f2) = a1.force2(a2)
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
    integrator.step()
  }
  
  def visualize(viz: Visualizer, disp: Atom => Visualizer.Sphere) {
    viz.setParticles(_atoms.map(disp(_)))
  }
}

