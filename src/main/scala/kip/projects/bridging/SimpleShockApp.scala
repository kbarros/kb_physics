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



object SimpleStress extends StressFn {
  val c = 0.5 // wave speed
  val rho0 = 0.5 // density
  
  def stress(defgrad: Double, energy: Double) = rho0*c*c*(defgrad - 1)
  
  def zeroTempEnergy(defgrad: Double) = {
    rho0*c*c*(defgrad-1)*(defgrad-1)/2
  }
}


class SimpleShockSolver(val L: Int, dx: Double, dt: Double) {
  // initial deformation gradient
  def defgrad0 = Array.tabulate(L)(i => if (i < L/2) 1.01 else 1.0)
  
  var sve1 = new ElastodynamicLaws1d(L, SimpleStress.rho0, defgrad0, SimpleStress)  
  var sve2 = new ElastodynamicLaws1d(L, SimpleStress.rho0, defgrad0, SimpleStress)  
  var time = 0.0
  
  def step() {
    CLIntegratorFirstOrder.iter(sve1, dx, dt)
    CLIntegratorSecondOrder.iter(sve2, dx, dt)
    time += dt
  }
  
  def mod(x: Int, y: Int): Int = {
    val rem = x % y
    if (rem < 0) rem + y else rem
  }
  
  def linearDefgradSolution(): Array[Double] = {
    val del_i = ((time * SimpleStress.c) / dx).toInt
    Array.tabulate(L) { i =>
      val i1 = mod(i+del_i, L)
      val i2 = mod(i-del_i, L)
      0.5*(defgrad0(i1) + defgrad0(i2))
    } 
  }
  
  def linearVelocitySolution() = {
    val del_i = ((time * SimpleStress.c) / dx).toInt
    Array.tabulate(L) { i =>
      val i1 = mod(i+del_i, L)
      val i2 = mod(i-del_i, L)
      SimpleStress.c * 0.5*(defgrad0(i1) - defgrad0(i2))
    } 
  }
  
  def linearEnergySolution(): Array[Double] = {
    val defgrad = linearDefgradSolution()
    val vel     = linearVelocitySolution()
    val rho0    = SimpleStress.rho0
    val ret = new Array[Double](L)
    for (i <- ret.indices) {
      ret(i) = SimpleStress.zeroTempEnergy(defgrad(i)) + rho0*vel(i)*vel(i)/2
    }
    ret
  }
}


object SimpleShockApp extends App {
  new Control(new SimpleShockApp(), "Shock wave")  
}

class SimpleShockApp extends Simulation {
  val defgradPlot = new Plot("Deformation Gradient")
  val energyPlot = new Plot("Energy")
  
  var sim: SimpleShockSolver = _
  
  def load(c: Control) {
    c.frameTogether("Plots", defgradPlot, energyPlot)

    params.add("L", 400)
    params.add("dx", 1.0)
    params.add("dt", 0.5)
    params.add("net energy")
  }
  
  def animate() {
    defgradPlot.registerLines("defgrad1", new PointSet(0, 1, sim.sve1.defgrad), Color.GRAY)
    defgradPlot.registerLines("defgrad2", new PointSet(0, 1, sim.sve2.defgrad), Color.BLACK)
    defgradPlot.registerLines("defgrad exact", new PointSet(0, 1, sim.linearDefgradSolution()), Color.RED)
    
    energyPlot.registerLines("energy 1", new PointSet(0, 1, sim.sve1.energy), Color.GRAY)
    energyPlot.registerLines("energy 2", new PointSet(0, 1, sim.sve2.energy), Color.BLACK)
    energyPlot.registerLines("energy exact", new PointSet(0, 1, sim.linearEnergySolution()), Color.RED)
    
    params.set("net energy", sim.sve2.energy.sum)
  }
  
  def clear() {
  }
  
  def run() {
    val L = params.iget("L")
    val dx = params.fget("dx")
    val dt = params.fget("dt")
    sim = new SimpleShockSolver(L=L, dx=dx, dt=dt)
    
    while (true) {
     Job.animate()
     for (i <- 0 until 1)
       sim.step()
     Thread.sleep(10)
    }
  }
}
