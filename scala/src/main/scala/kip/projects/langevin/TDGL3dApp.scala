package kip.projects.langevin

import scikit.util.Utilities.format
import scikit.graphics.dim3.Grid3D
import scikit.jobs.Control
import scikit.jobs.Job
import scikit.jobs.Simulation


object TDGL3dApp extends App {
  new Control(new TDGL3dApp(), "Allen Cahn Coarsening")
}

class TDGL3dApp extends Simulation {
  val grid = new Grid3D("Grid")
  var sim: TDGL2d = _
  
  def load(c: Control) {
    c.frame(grid)
    params.addm("dt", 5.0)
    params.addm("r", 1.0)
    params.add("L", 15.0)
    params.add("Random seed", 41)
    params.add("dx", 0.5)
    params.add("Time")
    params.add("Energy")
  }
  
  def animate() {
    sim.readParams(params)
    
    grid.setScale(-1, 1)
    val lp = sim.lp
    grid.registerData(lp, lp, lp, sim.phi)
    
    params.set("Time", format(sim.t))
    params.set("Energy", format(sim.freeEnergy))
  }
  
  def clear() {
    grid.clear()
  }
  
  def run() {
    sim = new TDGL2d(params, dimensions=3)
    Job.animate()
    
    while (true) {
      for (i <- 0 until 5)
	sim.simulate()
      Job.animate()
    }
  }
}

