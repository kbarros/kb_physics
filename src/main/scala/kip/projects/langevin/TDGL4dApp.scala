package kip.projects.langevin

import scikit.util.Utilities.format
import scikit.graphics.dim3.Grid3D
import scikit.jobs.Control
import scikit.jobs.Job
import scikit.jobs.Simulation

object TDGL4dApp extends App {
  new Control(new TDGL4dApp(), "Allen Cahn Coarsening")
}

class TDGL4dApp extends Simulation {
  val grids = Seq.tabulate(5)(i => new Grid3D("Grid"+i))
  var sim: TDGL = _

  def load(c: Control) {
    grids.foreach(g => c.frame(g))
    params.addm("dt", 5.0)
    params.addm("r", 1.0)
    params.add("L", 25.0)
    params.add("Random seed", 41)
    params.add("dx", 1.0)
    params.add("Time")
    params.add("Energy")
  }

  def animate() {
    sim.readParams(params)

    grids.foreach(_.setScale(-1, 1))
    val lp = sim.lp

    for ((g, i) <- grids.zipWithIndex) {
      val offset = (lp*lp*lp) * (lp*i/grids.size) 
      val sliceData = sim.phi.slice(offset, offset + (lp * lp * lp))
      g.registerData(lp, lp, lp, sliceData)
    }

    params.set("Time", format(sim.t))
    params.set("Energy", format(sim.freeEnergy))
  }

  def clear() {
    grids.foreach(_.clear())
  }

  def run() {
    sim = new TDGL(params, dimensions = 4)
    sim.randomize()
    Job.animate()

    while (true) {
      for (i <- 0 until 5)
        sim.simulate()
      Job.animate()
    }
  }
}

