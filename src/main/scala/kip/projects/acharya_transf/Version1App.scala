package kip.projects.acharya_transf

import scikit.jobs.Control
import scikit.jobs.Simulation
import scikit.graphics.dim2.Plot
import scikit.jobs.Job
import scikit.dataset.PointSet
import java.awt.Color

import math._

class Version1Sim(L: Int, dx: Double, dt: Double) {
  val omega = 1

  val A     = new Array[Double](L)
  val Ax    = new Array[Double](L)
  
  val f     = new Array[Double](L)
  val fdot  = new Array[Double](L)
  val fx    = new Array[Double](L)
  val fxdot = new Array[Double](L)
  
  val Tx = new Array[Double](L)
  
  { // initialize fields
    for (i <- 0 until L) {
      f(i) = 0.1*random 
      A(i) = f(i) //  0.1*random
    }
  }
  
  // T-tilde, derived from double well potential 
  def stress(fx: Double) = -fx + fx*fx*fx 
  
  def p(i: Int) = if (i == 0) L-1 else i-1
  def n(i: Int) = if (i == L-1) 0 else i+1
  
  def step() {
    // calculate force on f
    
    for (i <- 0 until L) {
      fx(i) = (f(n(i)) - f(i)) / dx
      fxdot(i) = (fdot(n(i)) - fdot(i)) / dx
      Ax(i) = (A(n(i)) - A(p(i))) / (2*dx)
    }
    
    val fxxx = Array.tabulate[Double](L) { i =>
      (fx(n(i)) - 2 * fx(i) + fx(p(i))) / (dx*dx) 
    }
    
    val Axx = Array.tabulate[Double](L) { i =>
      (A(n(i)) - 2 * A(i) + A(p(i))) / (dx*dx) 
    }
    
    val Ttilde = Array.tabulate[Double](L) { i =>
      -fx(i) + fx(i)*fx(i)*fx(i)
    }
    
    val effectiveStress = Array.tabulate[Double](L) { i =>
      val balance = 0.2
      Ttilde(i) - omega * ((1-balance) * Axx(i) + balance * fxxx(i)) 
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
      val m = 1
      val k = 500
      val constraint = 2*k*m*pow(abs(fx(i) - A(i)), 2*m-2) * (fx(i) - A(i))
      
      A(i) = 0.5 *(A(p(i)) + A(n(i))) + dt * abs(Ax(i)) * (constraint + omega*Axx(i)) + dt * fxdot(i)
    }
    
  }
}

object Version1App extends App {
  new Control(new Version1App(), "Phase transformation")  
}

class Version1App extends Simulation {
  val defgradPlot = new Plot("Deformation Gradient")
  val energyPlot = new Plot("Energy")
  
  var sim: Version1Sim = _
  
  def load(c: Control) {
    c.frameTogether("Plots", defgradPlot, energyPlot)

    params.add("L", 50)
    params.add("dx", 0.5)
    params.add("dt", 0.001)
    params.add("net energy")
  }
  
  def animate() {
    defgradPlot.registerLines("f", new PointSet(0, 1, sim.f), Color.BLUE)
    defgradPlot.registerLines("fx", new PointSet(0, 1, sim.fx), Color.GRAY)
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
    sim = new Version1Sim(L=L, dx=dx, dt=dt)
    
    while (true) {
     Job.animate()
     for (i <- 0 until 100)
       sim.step()
     Thread.sleep(10)
    }
  }
}
