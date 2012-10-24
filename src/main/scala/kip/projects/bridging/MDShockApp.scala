package kip.projects.bridging

import java.awt.Color

import scikit.dataset.PointSet
import scikit.graphics.dim2.Plot
import scikit.jobs.Control
import scikit.jobs.Job
import scikit.jobs.Simulation


object MDShockApp extends App {
  new Control(new MDShockApp(), "Shock wave")  
}

class MDShockApp extends Simulation {
  val defgradPlot = new Plot("Deformation Gradient")
  val energyPlot = new Plot("Energy")
  
  var sim: MDShockSolver = _
  
  // MDStressFn.checkZeroTempConsistency(1.1, 0.0001)
  
  def load(c: Control) {
    c.frameTogether("Plots", defgradPlot, energyPlot)
    
    params.add("L", 400)
    params.add("dx", 1.0)
    params.add("dt", 0.02)
    params.add("net energy")
    params.add("time")
  }
  
  def animate() {
    defgradPlot.registerLines("defgrad1", new PointSet(0, 1, sim.s1.defgrad), Color.GRAY)
    defgradPlot.registerLines("defgrad2", new PointSet(0, 1, sim.s2.defgrad), Color.BLACK)
    
    val potentialExact = Array.tabulate(sim.L)(i => MDStressFn.zeroTempEnergy(sim.s2.defgrad(i)))
    val potential1 = new Array[Double](sim.L)
    val potential2 = new Array[Double](sim.L)
    for (i <- 0 until sim.L) {
      val rho0 = MDStressFn.rho0
      val p1 = sim.s1.momentum(i)
      val p2 = sim.s2.momentum(i)
      potential1(i) = sim.s1.energy(i) - p1*p1/(2*rho0)
      potential2(i) = sim.s2.energy(i) - p2*p2/(2*rho0)
    }
    energyPlot.registerLines("energy 1", new PointSet(0, 1, potential1), Color.GRAY)
    energyPlot.registerLines("energy 2", new PointSet(0, 1, potential2), Color.BLACK)
    energyPlot.registerLines("energy exact", new PointSet(0, 1, potentialExact), Color.RED)
    
    params.set("net energy", sim.s2.energy.sum)
    params.set("time", sim.time)
  }
  
  def clear() {
  }
  
  def run() {
    val L = params.iget("L")
    val dx = params.fget("dx")
    val dt = params.fget("dt")
    sim = new MDShockSolver(L=L, dx=dx, dt=dt, defgradAmplitude=0.04, defgradWidth=5 /*5*dx*/)
    
    while (true) {
     Job.animate()
     for (i <- 0 until 10) {
       sim.step()
     }
     Thread.sleep(10)
    }
  }
}
