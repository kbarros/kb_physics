package kip.projects.acharya_transf

import scikit.jobs.Control
import scikit.jobs.Simulation
import scikit.graphics.dim2.Plot
import scikit.jobs.Job
import scikit.dataset.PointSet
import java.awt.Color

import math._


// Integrate.
//
//   phi_t + f = 0
//   f = - (phi_x)^2 (g(phi) + omega phi_xx)
//
// We transform this to a pseudo-conservation law,
//
//   u_t + f_x = 0
//   f = - u^2 (g(phi) + omega u_x) 
//   u = phi_x
//
// which can be integrated for u using Nessyahu-Tadmor. Then we need to reconstruct phi from
// u = phi_x.
// We use the approximate identity
//
//   \int dx phi(t+dt) = \int dx phi(t) - dt \int dx f(t+dt/2)
//
// To fix the constant C in:
//    phi = \int u dx + C 
//

class TadmorTestSim(val L: Int, val dx: Double, val dt: Double, val omega: Double) {
  val r = dt / dx

  var time = 0.0
  
  val phi0  = new Array[Double](L)
  val phi   = new Array[Double](L)
  val u     = new Array[Double](L)
  
  def m(i: Int) = (i-1+L) % L
  def p(i: Int) = (i+1) % L

  { // initialize fields
    for (i <- 0 until L) {
      val x = dx * i
      phi0(i) = 0.1 * cos(Pi * x / (dx*L/4))
      phi(i) = phi0(i)
    }
    
    for (i <- 0 until L) {
      u(i) = (phi(p(i)) - phi(m(i))) / (2*dx)
    }
  }
  
  def fluxFn(phi: Array[Double], u: Array[Double]) = Array.tabulate[Double](L) { i =>
    val ux = (u(p(i)) - u(m(i))) / (2*dx)
    - u(i) * u(i) * (phi(i) + omega * ux)
  }
  
  def nessyahuTadmorStep() {
    def minmod(a: Double, b: Double, c: Double): Double = {
      if (a > 0 && b > 0 && c > 0)
        min(a, min(b, c))
      else if (a < 0 && b < 0 && c < 0)
        max(a, max(b, c))
      else
        0
    }
    
    def deltaTilde(a: Array[Double], i: Int): Double = {
      val theta = 2
      minmod(
        theta * (a(p(i)) - a(i)),
        0.5 * (a(p(i)) - a(m(i))),
        theta * (a(i) - a(m(i)))
      )
    }
    
    def halfStep() {
      // f_i^n
      val f0 = fluxFn(phi, u)

      // u_i^{n+1/2}
      val u1 = Array.tabulate[Double](L) { i =>
        u(i) - (r/2) * deltaTilde(f0, i)
      }

      // f_i^{n+1/2}
      val f1 = fluxFn(phi, u1)

      // u_{i+1/2}^{n+1}
      val u2 = Array.tabulate[Double](L) { i =>
        (1/2f) * (u(i) + u(p(i))) +
        (1/8f) * (deltaTilde(u, i) - deltaTilde(u, p(i))) - r * (f1(p(i)) - f1(i))
      }
      
      for (i <- 0 until L) u(i) = u2(i)
      
      // build new phi
      val netPhi0 = phi.sum
      phi(0) = 0
      for (i <- 1 until L) {
        phi(i) = phi(i-1) + (dx/2) * (u(i-1) + u(i))
      }
      val netPhi1 = phi.sum
      
      val currentDelta = netPhi1 - netPhi0
      val exactDelta = - dt * f1.sum
      
      for (i <- 0 until L) {
        phi(i) += (exactDelta - currentDelta) / L
      }
    }
    
    // assign u_i = u_{i-1}
    def rotate(a: Array[Double]) {
      val an = a(L-1)
      for (i <- L-1 to 1 by -1)
        a(i) = a(m(i))
      a(0) = an
    }
    
    halfStep()
    halfStep()
    
    rotate(u)
    rotate(phi)
    
    time += 2*dt
  }
  
  def step() {
    nessyahuTadmorStep()
  }
}

object TadmorTestApp extends App {
  new Control(new TadmorTestApp(), "Phase transformation")  
}

class TadmorTestApp extends Simulation {
  val uPlot = new Plot("u")
  val phiPlot = new Plot("phi")
  
  var sim: TadmorTestSim = _
  
  def load(c: Control) {
    c.frameTogether("Plots", uPlot, phiPlot)

    params.add("L", 10)
    params.add("dx", 0.03)
    params.add("dt", 0.001)
    params.add("omega", 0.1)
    
    params.add("time")
  }
  
  def animate() {
    uPlot.registerLines("u", new PointSet(0, sim.dx, sim.u), Color.BLUE)
    phiPlot.registerLines("phi", new PointSet(0, sim.dx, sim.phi), Color.RED)
    
    params.set("time", sim.time)
  }
  
  def clear() {
    uPlot.clear()
    phiPlot.clear()
  }
  
  def run() {
    val L = params.fget("L")
    val dx = params.fget("dx")
    val Lp = (L/dx).toInt
    val dt = params.fget("dt")
    val omega = params.fget("omega")
    sim = new TadmorTestSim(L=Lp, dx=dx, dt=dt, omega=omega)
    
    while (true) {
     Job.animate()
     
     val lastTime = sim.time
     while (sim.time-lastTime < 0.049999) {
       sim.step()
     }
          
     Thread.sleep(10)
    }
  }
}
