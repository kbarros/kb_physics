//
// integrate coupled equations:
//   A_dot = |A_x| (1 - A_xx)
//   f_dot_dot = ...
//

package kip.projects.acharya_transf

import scikit.jobs.Control
import scikit.jobs.Simulation
import scikit.graphics.dim2.Plot
import scikit.jobs.Job
import scikit.dataset.PointSet
import java.awt.Color

import math._


class Version4Sim(val L: Int, val dx: Double, val dt: Double, val omega: Double, val v: Double) {
  var time = 0.0
  
  val A0 = 1   // preferred A0
  val f_viscosity = 1.0 // +1.0
  val f_drag = 0.01

  val A = new Array[Double](L)
  
  val f     = new Array[Double](L)
  val fdot  = new Array[Double](L)

  var f0 = 0.0
  var f1 = 0.0
  
  initD()
  fixFxToA()
  
  def fixFxToA() {
    f(0) = f0
    for (i <- 1 until L) {
      f(i) = f(i-1) + dx * A(i-1)
    }
    f1 = f(L-1)
  }
  
  def initD() {
    for (i <- 0 until L) {
      val x = dx * i
      A(i) = cos(Pi * x / (dx*L/4))
    }
  }
  
  def m(i: Int) = if (i == 0) 0 else i-1
  def p(i: Int) = if (i == L-1) L-1 else i+1
  
  
  def step() {
    updateA()
    updateF()
    
    // f0 = 0.1 * time // sin(1.0*time)
    time += dt
  }
  
  def fx(i: Int) = (f(p(i)) - f(i)) / dx
  
  def stress(i: Int) = { // d psi / d e
    val e = fx(i)-A(i)
    e * (e*e - A0*A0)
//    e
  }
  
  def tau(i: Int) = { // d eta / d A
    A(i) * (A(i)*A(i) - A0*A0)
  }
  
  def updateF() {
    val fxxx = Array.tabulate[Double](L) { i =>
      (fx(p(i)) - 2 * fx(i) + fx(m(i))) / (dx*dx) 
    }
    
    val effectiveStress = Array.tabulate[Double](L) { i =>
      stress(i) + f_viscosity * fx(i) // - f_viscosity * fxxx(i) 
    }
    
    // velocity-verlet to update f and fdot 
    for (i <- 0 until L) {
      val force = (effectiveStress(i) - effectiveStress(m(i))) / dx 
      fdot(i) += dt * (force - f_drag * fdot(i))
      f(i) += dt * fdot(i)
    }
    f(0) = f0
    f(L-1) = f1
  }
  
  def updateA() {
    def sqr(x: Double) = x*x
    
    def force(i: Int) = stress(i) - tau(i) // Axx term handled implicitly
    
    val absAx = Array.tabulate[Double](L) { i =>
      val Axp = (A(p(i)) - A(i)) / dx
      val Axm = (A(i) - A(m(i))) / dx
      val Axx = (A(p(i)) - 2 * A(i) + A(m(i))) / (dx*dx) 

      // A_t = |A_x| F 
      val F = force(i) + omega * Axx
      
      if (F < 0) {
        sqrt(sqr(max(Axm, 0)) +  sqr(min(Axp, 0)))
      }
      else {
        sqrt(sqr(max(Axp, 0)) + sqr(min(Axm, 0)))
      }
      
      1
    }
    
    import no.uib.cipr.matrix._

    val rhs = new DenseVector(Array.tabulate[Double](L) { i =>
      // can't do periodic boundaries (would be non tridiagonal) 
      val extra = if (i == 0 || i == L-1) omega * ((dt * absAx(i)) / (dx*dx)) * A(i) else 0
      A(i) + dt * force(i) * absAx(i) + extra
    })

    val tdmat = new TridiagMatrix(L)
    val ld = tdmat.getSubDiagonal()
    val d = tdmat.getDiagonal()
    val ud = tdmat.getSuperDiagonal()
    for (i <- 0 until L) {
      val pf = omega * dt * absAx(i) / (dx*dx)

      if (i > 0) ld(i-1) = - pf
      d(i)  = 1 + 2*pf
      if (i < L-1) ud(i) = - pf
    }
    
    val x = new DenseVector(A)
    tdmat.solve(rhs, x)
    
    for (i <- 1 until L-1) {
      A(i) = x.get(i)
    }
    
    // enforce boundary condition A_x = 0
    A(0) = A(1)
    A(L-1) = A(L-2)
  }
}

object Version4App extends App {
  new Control(new Version4App(), "Phase transformation")  
}

class Version4App extends Simulation {
  val defgradPlot = new Plot("Deformation Gradient")
  
  var sim: Version4Sim = _
  
  def load(c: Control) {
    c.frameTogether("Plots", defgradPlot)
    
    params.add("L", 12.0)
    params.add("dt", 0.001)
    params.add("dx", 0.01)
    params.add("omega", 0.01)
    params.add("v", 0.0)
    params.add("time")
  }
  
  def animate() {
    defgradPlot.registerLines("A", new PointSet(0, sim.dx, sim.A), Color.RED)
    defgradPlot.registerLines("f", new PointSet(0, sim.dx, sim.f), Color.BLUE)
    val fx = Array.tabulate(sim.L)(sim.fx(_))
    defgradPlot.registerLines("fx", new PointSet(0, sim.dx, fx), Color.PINK)
    
    params.set("time", sim.time)
  }
  
  def clear() {
    defgradPlot.clear()
  }
  
  def run() {
    val L = params.fget("L")
    val dx = params.fget("dx")
    val dt = params.fget("dt")
    val omega = params.fget("omega")
    val v = params.fget("v")
    
    val Lp = {
      val t = (L/dx).toInt
      t - t % 12
    }
    
    sim = new Version4Sim(L=Lp, dx=dx, dt=dt, omega=omega, v=v)
    
    while (true) {

     Job.animate()
     
     val lastTime = sim.time
     while (sim.time-lastTime < 0.1) {
       sim.step()
     }
     
     Thread.sleep(10)
    }
  }
}