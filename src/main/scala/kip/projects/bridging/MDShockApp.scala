package kip.projects.bridging

import scikit.dataset.PointSet
import scikit.jobs.Control
import scikit.jobs.Job
import scikit.jobs.Simulation
import scikit.jobs.params.DoubleValue
import scikit.jobs.params.Parameters
import scikit.numerics.Math2;
import scala.util.Random
import scala.math.sqrt

import java.awt.Color
import scikit.graphics.dim2.{Scene2D, Geom2D, Plot}



object MDStressFn extends StressFn {
  val md = new SubMD2d(ncols=2, nrows=2, a=1.12246, dt=0.01)
  val rho0 = md.density0
  
  def stress(defgrad: Double, energy: Double): Double = {
    md.applyDeformationGradient(defgrad)
    md.initializeCrystalPositions()
    val zeroTemperature = true
    if (zeroTemperature) {
      md.applyKineticEnergy(0)
      md.convertCauchyToFirstPiolaKirchoff(md.virialStress())
    }
    else {
      md.applyKineticEnergy(energy*md.volume0 - md.potentialEnergy())
      md.convertCauchyToFirstPiolaKirchoff(md.virialStress()) // TODO average data
    }
  }
  
  def zeroTempEnergy(defgrad: Double) = {
    md.applyDeformationGradient(defgrad)
    md.initializeCrystalPositions()
    md.applyKineticEnergy(0)
    md.potentialEnergy() / md.volume0
  }
}


class MDShockSolver(val L: Int, dx: Double, dt: Double, defgradAmplitude: Double, defgradWidth: Double) {
  // initial deformation gradient
//  def defgrad0 = Array.tabulate(L)(i => 1.0 + 0.04*(math.tanh(0.1*(i - L/2)) + 1.0))
  
  import math.tanh
  val (a, w) = (defgradAmplitude, defgradWidth)
  def defgrad0 = Array.tabulate(L)(i => 1.0 + (a/2)*(tanh((dx/w)*(i-L/4)) - tanh((dx/w)*(i-3*L/4))))
  
  var sve1 = new ElastodynamicLaws1d(L, MDStressFn.rho0, defgrad0, MDStressFn)  
  var sve2 = new ElastodynamicLaws1d(L, MDStressFn.rho0, defgrad0, MDStressFn)  
  var time = 0.0
  
  def step() {
    CLIntegratorFirstOrder.iter(sve1, dx, dt)
    CLIntegratorSecondOrder.iter(sve2, dx, dt)
    time += dt
  }
}


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
  }
  
  def animate() {
    defgradPlot.registerLines("defgrad1", new PointSet(0, 1, sim.sve1.defgrad), Color.GRAY)
    defgradPlot.registerLines("defgrad2", new PointSet(0, 1, sim.sve2.defgrad), Color.BLACK)
    
    val potentialExact = Array.tabulate(sim.L)(i => MDStressFn.zeroTempEnergy(sim.sve2.defgrad(i)))
    val potential1 = new Array[Double](sim.L)
    val potential2 = new Array[Double](sim.L)
    for (i <- 0 until sim.L) {
      val rho0 = MDStressFn.rho0
      val p1 = sim.sve1.momentum(i)
      val p2 = sim.sve2.momentum(i)
      potential1(i) = sim.sve1.energy(i) - p1*p1/(2*rho0)
      potential2(i) = sim.sve2.energy(i) - p2*p2/(2*rho0)
    }
    energyPlot.registerLines("energy 1", new PointSet(0, 1, potential1), Color.GRAY)
    energyPlot.registerLines("energy 2", new PointSet(0, 1, potential2), Color.BLACK)
    energyPlot.registerLines("energy exact", new PointSet(0, 1, potentialExact), Color.RED)
    
    params.set("net energy", sim.sve2.energy.sum)
  }
  
  def clear() {
    println("clear")
  }
  
  def run() {
    val L = params.iget("L")
    val dx = params.fget("dx")
    val dt = params.fget("dt")
    sim = new MDShockSolver(L=L, dx=dx, dt=dt, defgradAmplitude=0.04, defgradWidth=5*dx)
    
    while (true) {
     Job.animate()
     for (i <- 0 until 1)
       sim.step()
     Thread.sleep(10)
    }
  }
}
