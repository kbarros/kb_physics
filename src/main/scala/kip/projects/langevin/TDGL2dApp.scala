package kip.projects.langevin

import scikit.util.Utilities.format
import scikit.graphics.dim2.Grid
import scikit.jobs.Control
import scikit.jobs.Job
import scikit.jobs.Simulation


object TDGL2dApp extends App {
  new Control(new TDGL2dApp(), "Allen Cahn Coarsening")
}

class TDGL2dApp extends Simulation {
  val grid = new Grid("Grid")
  var sim: TDGL = _

  def load(c: Control) {
    c.frame(grid)
    params.addm("dt", 1.0)
    params.addm("r", 1.0);
    params.add("L", 500.0)
    params.add("Random seed", 0)
    params.add("dx", 1.0)
    params.add("Time")
    params.add("Energy")
  }
  
  def animate() {
    sim.readParams(params)
    
    grid.setScale(-1, 1)
    val lp = sim.lp
    grid.registerData(lp, lp, sim.phi)

    params.set("Time", format(sim.t))
    params.set("Energy", format(sim.freeEnergy))
  }
  
  def clear() {
    grid.clear()
  }
  
  def run() {
    sim = new TDGL(params, dimensions=2)
    sim.randomize()
    Job.animate()
    
    while (true) {
      for (i <- 0 until 10)
	sim.simulate()
      Job.animate()
    }
  }
}
