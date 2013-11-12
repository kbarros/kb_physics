package kip.projects.acharya_transf

import scikit.jobs.Control
import scikit.jobs.Simulation
import scikit.graphics.dim2.Plot
import scikit.jobs.Job
import scikit.dataset.PointSet
import java.awt.Color
import math._
import kip.math.Math
import Math._



// Integrate.
//
//   phi1_t + f1 = 0
//   phi2_t + f2 = 0
//
//   f1 = - S phi1_x
//   f2 = - S phi2_x
//   S = (phi1 + eps phi1_xx) phi1_x + (phi2 + eps phi2_xx) phi2_x
//
// We transform this to a pseudo-conservation law,
//
//   u1_t + f1_x = 0
//   u2_t + f2_x = 0
//
//   f1 = - S u1
//   f2 = - S u2
//   S = (phi1 + eps u1_x) u1 + (phi2 + eps u2_x) u2
//
// which can be integrated for u using Nessyahu-Tadmor. Then we need to
// reconstruct phi1 (and phi2) from u1 = phi1_x
//
// We use the approximate identity
//
//   \int dx phi1(t+dt) = \int dx phi1(t) - dt \int dx f1(t+dt/2)
//
// To fix the constant C in:
//    phi1 = \int u1 dx + C 
//

class TartarTestSim(val L: Int, val dx: Double, val dt: Double, val omega: Double) {
  val r = dt / dx

  var time = 0.0
  
  val phi1   = new Array[Double](L)
  val phi2   = new Array[Double](L)
  val u1     = new Array[Double](L)
  val u2     = new Array[Double](L)
  
  def m(i: Int) = (i-1+L) % L
  def p(i: Int) = (i+1) % L

  // initialize fields
  def initA() { 
    for (i <- 0 until L) {
      val x = dx * i
      phi1(i) = 1 + 0.1 * sin(Pi * x / (dx*L/8))// + 0.005 * sin(Pi*x / (dx*L/64))
      phi2(i) = 1 + 0.1 * cos(Pi * x / (dx*L/8))
//      phi1(i) = 0.1 * exp(- sqr(x - L*dx/2) / sqr(0.5))
//      phi2(i) = 0.05 * exp(- sqr(x - L*dx/2) / sqr(2)) * cos(Pi * x / (dx*L/17))
    }
  }
  
  def initB() {
    for (i <- 0 until L) {
      val x = dx * i
      
      val a = 0
      val b = 1
      phi1(i) = i*4 / L match {
      case 0 => a
      case 1 => {
        val theta = 4.0*i/L - 1.0
        b*theta + a*(1-theta)
      }
      case 2 => b
      case 3 => {
        val theta = 4.0*i/L - 3.0
        a*theta + b*(1-theta)
      }
      }
    }
    for (i <- 0 until L) {
      phi2(i) = phi1((i+L/8)%L)
    }
  }

  initB()
  for (i <- 0 until L) {
    u1(i) = (phi1(p(i)) - phi1(m(i))) / (2*dx)
    u2(i) = (phi2(p(i)) - phi2(m(i))) / (2*dx)
  }
  
  // force = - d energy / d phi
  def force(phi1: Double, phi2: Double) = (+phi1, +phi2)
  def energy(phi1: Double, phi2: Double) = 0.5 * (sqr(phi1) + sqr(phi2))
  
  
  def fluxFns(phi1: Array[Double],
              phi2: Array[Double],
              u1: Array[Double],
              u2: Array[Double]) = {
    val fs = Array.tabulate[(Double, Double)](L) { i =>
      val x = i * dx
      val h = 0 // 0.01 * exp(- sqr(x - L*dx/2) / sqr(0.5))
      
      val u1x = (u1(p(i)) - u1(m(i))) / (2*dx)
      val u2x = (u2(p(i)) - u2(m(i))) / (2*dx)
      
      //val g1 = -phi1(i) // phi1(i) - cube(phi1(i))
      //val g2 = -phi2(i)
      val (g1, g2) = force(phi1(i), phi2(i))
      
      val S = (g1 + omega * u1x) * u1(i) + (g2 + omega * u2x) * u2(i)
      (- S * u1(i) + h, - S * u2(i) + h)
    }
    val (f1, f2) = fs.unzip
    (f1.toArray, f2.toArray)
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
    
    def buildPhi(f: Array[Double], u: Array[Double], phi: Array[Double]) {
      val netPhi0 = phi.sum
      phi(0) = 0
      for (i <- 1 until L) {
        phi(i) = phi(i-1) + (dx/2) * (u(i-1) + u(i))
      }
      val netPhi1 = phi.sum
      
      val currentDelta = netPhi1 - netPhi0
      val exactDelta = - dt * f.sum
      
      for (i <- 0 until L) {
        phi(i) += (exactDelta - currentDelta) / L
      }
    }
    
    def halfStep() {
      // f_i^n
      val (fa1, fa2) = fluxFns(phi1, phi2, u1, u2)

      // u_i^{n+1/2}
      val ub1 = Array.tabulate[Double](L) { i =>
        u1(i) - (r/2) * deltaTilde(fa1, i)
      }
      val ub2 = Array.tabulate[Double](L) { i =>
        u2(i) - (r/2) * deltaTilde(fa2, i)
      }
      
      // f_i^{n+1/2}
      val (fb1, fb2) = fluxFns(phi1, phi2, ub1, ub2)
      
      // u_{i+1/2}^{n+1}
      val uc1 = Array.tabulate[Double](L) { i =>
        (1/2f) * (u1(i) + u1(p(i))) +
        (1/8f) * (deltaTilde(u1, i) - deltaTilde(u1, p(i))) - r * (fb1(p(i)) - fb1(i))
      }
      val uc2 = Array.tabulate[Double](L) { i =>
        (1/2f) * (u2(i) + u2(p(i))) +
        (1/8f) * (deltaTilde(u2, i) - deltaTilde(u2, p(i))) - r * (fb2(p(i)) - fb2(i))
      }
      
      for (i <- 0 until L) {
        u1(i) = uc1(i)
        u2(i) = uc2(i)
      }
      
      buildPhi(fb1, u1, phi1)
      buildPhi(fb2, u2, phi2)
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
    
    rotate(u1)
    rotate(u2)
    rotate(phi1)
    rotate(phi2)
    
    time += 2*dt
  }
  
  def step() {
    nessyahuTadmorStep()
  }
}


object TartarTestApp extends App {
  new Control(new TartarTestApp(), "Phase transformation")  
}

class TartarTestApp extends Simulation {
  val uPlot = new Plot("u")
  val phiPlot = new Plot("phi")
  
  var sim: TartarTestSim = _
  
  def load(c: Control) {
    c.frameTogether("Plots", uPlot, phiPlot)

    params.add("L", 10)
    params.add("dx", 0.01)
    params.add("dt", 0.001)
    params.add("omega", 0.01)
    
    params.add("time")
  }
  
  def animate() {
    uPlot.registerLines("u1", new PointSet(0, sim.dx, sim.u1), Color.BLUE)
    uPlot.registerLines("u2", new PointSet(0, sim.dx, sim.u2), Color.RED)
    phiPlot.registerLines("phi1", new PointSet(0, sim.dx, sim.phi1), Color.BLUE)
    phiPlot.registerLines("phi2", new PointSet(0, sim.dx, sim.phi2), Color.RED)
    
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
    sim = new TartarTestSim(L=Lp, dx=dx, dt=dt, omega=omega)
    
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
