//
// integrate A_dot = |A_x| (1 - A_xx)
//
// with pointy initial conditions
//

package kip.projects.acharya_transf

import scikit.jobs.Control
import scikit.jobs.Simulation
import scikit.graphics.dim2.Plot
import scikit.jobs.Job
import scikit.dataset.PointSet
import java.awt.Color

import math._


class Version3Sim(val L: Int, val dx: Double, val dt: Double) {
  var time = 0.0
  
  val A0 = new Array[Double](L)
  val A = new Array[Double](L)
  
  initB()
  
  def initB() {
    for (i <- 0 until L) {
      val Lp = L / 3
      require(L % 3 == 0)
      
      val x = (i / Lp) match {
        case 0 | 2 => {
          dx*Lp/2
        }
        case 1 => {
          val ip = (i % Lp) - Lp/2
          ip*dx
        }
      }
      A0(i) = x*x/2
      A(i) = A0(i)
    }

  }
  
  def initA() { // initialize spiky fields
    for (i <- 0 until L) {
      val nb = 4
      val Lp = L / nb
      
      val flat = false // (i / Lp) % 2 == 0
      val ip = (if (flat) 0 else (i % Lp)) - Lp/2
      
      val x = ip*dx
      
      A0(i) = x*x/2 / 2
      A(i) = A0(i)
    }
  }
  
  def p(i: Int) = if (i == 0) L-1 else i-1
  def n(i: Int) = if (i == L-1) 0 else i+1
  
  def laxWendroffStep() {
    
    // Ap(i) = A^{n}(i+1/2)
    val Ap = Array.tabulate[Double](L) { i =>
      (A(n(i)) + A(i)) / 2
    }
    
    // Axp(i) == Ax^{n}(i+1/2)
    val Axp = Array.tabulate[Double](L) { i =>
      (A(n(i)) - A(i)) / dx
    }
    
    // Axxp(i) == Axx^{n}(i+1/2)
    val Axxp = Array.tabulate[Double](L) { i =>
      (A(n(n(i))) - A(n(i)) - A(i) + A(p(i))) / (2*dx*dx) 
    }
    
    // Bp(i) = A^{n+1/2}(i+1/2)
    val Bp = Array.tabulate[Double](L) { i =>
      Ap(i) + (dt/2) * abs(Axp(i)) * (1 - Axxp(i)) 
    }
    
    // Bx(i) = Ax^{n+1/2}(i)
    val Bx = Array.tabulate[Double](L) { i =>
      (Bp(i) - Bp(p(i))) / dx
    }

    // Bxx(i) = Axx^{n+1/2}(i)
    val Bxx = Array.tabulate[Double](L) { i =>
      (Bp(n(i)) - Bp(i) - Bp(p(i)) + Bp(p(p(i)))) / dx
    }
    
    // Lax-Wendroff step to update A
    for (i <- 0 until L) {
      A(i) = A(i) + dt * abs(Bx(i)) * (1 - Bxx(i))
    }
    
    time += dt
  }

  def laxFriedrichsStep() {
    val Ax = Array.tabulate[Double](L) { i =>
      (A(n(i)) - A(p(i))) / (2*dx)
    }
    
    val Axx = Array.tabulate[Double](L) { i =>
      (A(n(i)) - 2 * A(i) + A(p(i))) / (dx*dx) 
    }
    
    for (i <- 0 until L) {
      A(i) = 0.5*(A(p(i)) + A(n(i))) + dt * abs(Ax(i)) * (1 - Axx(i))
    }
    
    time += dt
  }

}

object Version3App extends App {
  new Control(new Version3App(), "Phase transformation")  
}

class Version3App extends Simulation {
  val defgradPlot = new Plot("Deformation Gradient")
  val energyPlot = new Plot("Energy")
  
  var sim: Version3Sim = _
  var sim2: Version3Sim = _
  
  def load(c: Control) {
    c.frameTogether("Plots", defgradPlot, energyPlot)
    
    params.add("L", 10.0)
    params.add("dt = dx * ", 0.1)
    params.add("dx", 0.2)
    params.add("net energy")
  }
  
  def animate() {
    defgradPlot.registerLines("A0", new PointSet(0, sim.dx, sim.A0), Color.BLUE)
    defgradPlot.registerLines("A1", new PointSet(0, sim.dx, sim.A), Color.RED)
    defgradPlot.registerLines("A2", new PointSet(0, sim2.dx, sim2.A), Color.PINK)
  }
  
  def clear() {
    defgradPlot.clear()
  }
  
  def run() {
    val L = params.fget("L")
    val dx = params.fget("dx")
    
    val Lp = {
      val t = (L/dx).toInt
      t - t % 12
    }
    
    val dt = params.fget("dt = dx * ") * dx
    
    printf("dx = %f, dt = %f\n", dx, dt)
    
    sim = new Version3Sim(L=Lp, dx=dx, dt=dt)
    sim2 = new Version3Sim(L=Lp/2, dx=dx*2, dt=dt)
    
    while (true) {

     Job.animate()
     
     val lastTime = sim.time
     while (sim.time-lastTime < 0.1) {
//       sim.laxWendroffStep()
       sim.laxFriedrichsStep()
       sim2.laxFriedrichsStep()
     }
     
     Thread.sleep(10)
    }
  }
}