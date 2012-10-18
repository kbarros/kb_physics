package kip.projects.acharya_transf

import scikit.jobs.Control
import scikit.jobs.Simulation
import scikit.graphics.dim2.Plot
import scikit.jobs.Job
import scikit.dataset.PointSet
import java.awt.Color

import math._

class Version2Sim(L: Int, dx: Double, dt: Double) {
  var time = 0.0
  
  val omega = 1
  val f_viscosity = 0.01
  
  val f     = new Array[Double](L)
  val fx    = new Array[Double](L)
  
  val A     = new Array[Double](L)
  
  val fdot  = new Array[Double](L)
  
  { // initialize fields
    for (i <- 0 until L) {
      // val lambda = L*dx / 4
      // f(i) = sin(2*Pi * i*dx/lambda)

      val nb = 4
      val Lp = L / nb
      
      val flat = false // (i / Lp) % 2 == 0
      val ip = (if (flat) 0 else (i % Lp)) - Lp/2
      
      val x = ip*dx
      
      f(i) = x*x/2
      A(i) = f(i)
    }
  }
  
  // T-tilde, derived from double well potential
  
  def stress(e: Double) = -e + e*e*e
  
  def p(i: Int) = if (i == 0) L-1 else i-1
  def n(i: Int) = if (i == L-1) 0 else i+1
  
  def step() {
    // calculate force on f
    
    for (i <- 0 until L) {
      fx(i) = (f(n(i)) - f(i)) / dx
    }
    
    val Ax = Array.tabulate[Double](L) { i =>
      (A(n(i)) - A(p(i))) / (2*dx)
    }
    
    val fxxx = Array.tabulate[Double](L) { i =>
      (fx(n(i)) - 2 * fx(i) + fx(p(i))) / (dx*dx) 
    }
    
    val Axx = Array.tabulate[Double](L) { i =>
      (A(n(i)) - 2 * A(i) + A(p(i))) / (dx*dx) 
    }
    
    val Ttilde = Array.tabulate[Double](L) { i =>
      stress(fx(i) - A(i))
    }
    
    val effectiveStress = Array.tabulate[Double](L) { i =>
      Ttilde(i) - f_viscosity * fxxx(i) 
    }
    
    // velocity-verlet to update f and fdot 
    for (i <- 0 until L) {
      val drag = 0.01
      val force = (effectiveStress(i) - effectiveStress(p(i))) / dx 
      fdot(i) += dt * (force - drag * fdot(i))
      f(i) += dt * fdot(i)
    }
    
    // Lax-Friedrichs step to update A
    for (i <- 0 until L) {
      A(i) = 0.5 *(A(p(i)) + A(n(i))) + dt * abs(Ax(i)) * (-1 /* Ttilde(i) */ + omega*Axx(i))
    }
    
    time += dt
  }
}

object Version2App extends App {
  new Control(new Version2App(), "Phase transformation")  
}

class Version2App extends Simulation {
  val defgradPlot = new Plot("Deformation Gradient")
  val energyPlot = new Plot("Energy")
  
  var sim: Version2Sim = _
  
  def load(c: Control) {
    c.frameTogether("Plots", defgradPlot, energyPlot)

    params.add("L", 5e3.toInt)
    params.add("dx", 5e-3)
    params.add("dt", 1e-6)
    params.add("net energy")
  }
  
  def animate() {
    defgradPlot.registerLines("f", new PointSet(0, 1, sim.f), Color.BLUE)
//    defgradPlot.registerLines("fx", new PointSet(0, 1, sim.fx), Color.GRAY)
    defgradPlot.registerLines("A", new PointSet(0, 1, sim.A), Color.RED)
//    defgradPlot.registerLines("Ax", new PointSet(0, 1, sim.Ax), Color.PINK)

    // energyPlot.registerLines("df/dt", new PointSet(0, 1, sim.fdot), Color.BLACK)
    
    // params.set("net energy", sim.sve2.energy.sum)
  }
  
  def clear() {
  }
  
  def run() {
    val L = params.iget("L")
    val dx = params.fget("dx")
    val dt = params.fget("dt")
    sim = new Version2Sim(L=L, dx=dx, dt=dt)
    
    while (true) {
     Job.animate()
     
//     for (i <- 0 until 1)
     val lastTime = sim.time
     while (sim.time-lastTime < 0.001) {
       sim.step()
     }
     
     Thread.sleep(10)
    }
  }
}
