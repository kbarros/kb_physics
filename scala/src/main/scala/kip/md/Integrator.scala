package kip.md

import scala.math._

class World {
  var periodic = true
  var temperature = 1.0
  var integrator = new Verlet(this, 0.1)

  private var _atoms = Seq[Atom]()
  private var _globalCutoff: Double = _
  private var _Lx: Double = _
  private var _Ly: Double = _
  private var _grid: PointGrid2d[Atom] = _
  private var _interactionRange: Double = _
  
  def setSize(Lx: Double, Ly: Double) {
    _Lx = Lx
    _Ly = Ly
    val atomsPerCell = 4
    val cols = sqrt(_atoms.size / atomsPerCell).toInt
    _grid = new PointGrid2d(_Lx, cols, periodic)
  }
  
  def setAtoms(atoms: Seq[Atom]) {
    _atoms = atoms
    _interactionRange = {
      val is = atoms.toSet[Atom].flatMap(_.tag.interactions)
      println("Num interactions = "+is.size)
      val x = for (i <- is; j <- i.compatibleInteractions(is)) yield i.cutoff(j)
      0.0
    }
      
  }

  def potentialEnergy {
    var V = 0.0
    _grid.loadPoints(_atoms)
    for (a1 <- _atoms) {
      for (a2 <- _grid.pointOffsetsWithinRange(a1, _interactionRange)) {
	if (a1 != a2)
	  V += a1.potential(a2)/2. // divisor of 2 corrects for double counting
      }
      // accumulate accelerations due to external forces
      // V += a1.potential();
    }
    V
  }
	
	
	// public double kineticEnergy() {
	// 	double K = 0;
	// 	for (Pt p : particles) {
	// 		double M = p.tag.mass;
	// 		K += 0.5*M*(p.vx*p.vx+p.vy*p.vy);
	// 	}
	// 	return K;
	// }

}

trait Integrator {
  var dt: Double
  def step()
  def time: Double
}

class Verlet(world: World, var dt: Double) extends Integrator {
  var time = 0.0
  
  def step() {
    
  }
}
