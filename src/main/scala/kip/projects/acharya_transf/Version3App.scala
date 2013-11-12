//
// integrate A_t = |A_x| (-v + omega A_xx)
//
//

package kip.projects.acharya_transf

import scikit.jobs.Control
import scikit.jobs.Simulation
import scikit.graphics.dim2.Plot
import scikit.jobs.Job
import scikit.dataset.PointSet
import java.awt.Color

import math._


class Version3Sim(val L: Int, val dx: Double, val dt: Double, val omega: Double, val v: Double) {
  var time = 0.0
  
  val A0 = new Array[Double](L)
  val A = new Array[Double](L)
  
  initD()
  
  def initD() {
    for (i <- 0 until L) {
      val x = dx * i
      A0(i) = cos(Pi * x / (dx*L/4))
      A(i) = A0(i)
    }
  }
  
  def initC() {
    for (i <- 0 until L) {
      val Lp = L / 3
      require(L % 3 == 0)
      
      val x = (i / Lp) match {
        case 0 => 0
        case 1 => dx * (i % Lp)
        case 2 => dx * Lp
      }
      
      A0(i) = -x*x/2
      A(i) = A0(i)
    }
  }

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
      A0(i) = -x*x/2
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
  
  def m(i: Int) = if (i == 0) 0 else i-1
  def p(i: Int) = if (i == L-1) L-1 else i+1
  
  
  def laxFriedrichsStep() {
    val Ax = Array.tabulate[Double](L) { i =>
      (A(p(i)) - A(m(i))) / (2*dx)
    }
    
    val Axx = Array.tabulate[Double](L) { i =>
      (A(p(i)) - 2 * A(i) + A(m(i))) / (dx*dx) 
    }
    
    for (i <- 0 until L) {
      A(i) = 0.5*(A(m(i)) + A(p(i))) + dt * abs(Ax(i)) * (- 1 + Axx(i))
    }
    
    time += dt
  }
  
  def laxWendroffStep() {
    
    // Ap(i) = A^{n}(i+1/2)
    val Ap = Array.tabulate[Double](L) { i =>
      (A(p(i)) + A(i)) / 2
    }
    
    // Axp(i) == Ax^{n}(i+1/2)
    val Axp = Array.tabulate[Double](L) { i =>
      (A(p(i)) - A(i)) / dx
    }
    
    // Axxp(i) == Axx^{n}(i+1/2)
    val Axxp = Array.tabulate[Double](L) { i =>
      (A(p(p(i))) - A(p(i)) - A(i) + A(m(i))) / (2*dx*dx) 
    }
    
    // Bp(i) = A^{n+1/2}(i+1/2)
    val Bp = Array.tabulate[Double](L) { i =>
      Ap(i) + (dt/2) * abs(Axp(i)) * (- 1 + Axxp(i)) 
    }
    
    // Bx(i) = Ax^{n+1/2}(i)
    val Bx = Array.tabulate[Double](L) { i =>
      (Bp(i) - Bp(m(i))) / dx
    }

    // Bxx(i) = Axx^{n+1/2}(i)
    val Bxx = Array.tabulate[Double](L) { i =>
      (Bp(p(i)) - Bp(i) - Bp(m(i)) + Bp(m(m(i)))) / (2*dx*dx)
    }
    
    // Lax-Wendroff step to update A
    for (i <- 1 until L-1) {
      A(i) = A(i) + dt * abs(Bx(i)) * (- 1 + Bxx(i))
    }
    
    time += dt
  }

  def semiImplicitEulerStep() {
    
    val absAx = Array.tabulate[Double](L) { i =>
      if (false) {
        val Ax1 = (A(p(i)) - A(i)) / dx
        val Ax2 = (A(i) - A(m(i))) / dx
        if (Ax1 > 0 != Ax2 > 0)
          0
        else
          abs(0.5 * (Ax1 + Ax2))
      }
      
      else {          
        val Ax1 = (A(p(i)) - A(i)) / dx
        val Ax2 = (A(i) - A(m(i))) / dx
        0.5 * (abs(Ax1) + abs(Ax2))
      }
    }
        
    import no.uib.cipr.matrix._

    val rhs = new DenseVector(Array.tabulate[Double](L) { i =>
      // don't use implicit scheme for wrap-around data (otherwise, non tridiagonal) 
      val extra = if (i == 0 || i == L-1) omega * ((dt * absAx(i)) / (dx*dx)) * A(i) else 0
      A(i) - dt * v * (i*dx-6) * absAx(i) + extra
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
    
    time += dt
  }

  def upwindStep() {
    def sqr(x: Double) = x*x
    
    val Axx = Array.tabulate[Double](L) { i =>
      (A(p(i)) - 2 * A(i) + A(m(i))) / (dx*dx) 
    }

    val Axp = Array.tabulate[Double](L) { i =>
      (A(p(i)) - A(i)) / dx
    }

    val Axm = Array.tabulate[Double](L) { i =>
      (A(i) - A(m(i))) / dx
    }

    val deltaP = Array.tabulate[Double](L) { i =>
      sqrt(sqr(max(Axm(i), 0)) +  sqr(min(Axp(i), 0)))
    }
    
    val deltaM = Array.tabulate[Double](L) { i =>
      sqrt(sqr(max(Axp(i), 0)) + sqr(min(Axm(i), 0)))
    }

    val F = Array.tabulate[Double](L) { i =>
      v*(i*dx-6) - omega * Axx(i)
    }
    
    for (i <- 1 until L-1) {
      A(i) = A(i) - dt * (max(F(i), 0) * deltaP(i) + min(F(i), 0) * deltaM(i))
    }

    time += dt
  }
  
  
  def implicitUpwindStep() {
    def sqr(x: Double) = x*x
    
    def vel(i: Int) = 1 // v * (i*dx - 6)
    
    val absAx = Array.tabulate[Double](L) { i =>
      val Axp = (A(p(i)) - A(i)) / dx
      val Axm = (A(i) - A(m(i))) / dx
      val Axx = (A(p(i)) - 2 * A(i) + A(m(i))) / (dx*dx) 

      // A_t = - |A_x| F 
      val F = vel(i) - omega * Axx
      
      if (F > 0) {
        sqrt(sqr(max(Axm, 0)) +  sqr(min(Axp, 0)))
      }
      else {
        sqrt(sqr(max(Axp, 0)) + sqr(min(Axm, 0)))
      }
    }
    
    import no.uib.cipr.matrix._

    val rhs = new DenseVector(Array.tabulate[Double](L) { i =>
      // can't do periodic boundaries (would be non tridiagonal) 
      val extra = if (i == 0 || i == L-1) omega * ((dt * absAx(i)) / (dx*dx)) * A(i) else 0
      A(i) - dt * vel(i) * absAx(i) + extra
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
    
    time += dt
  }
}

object Version3App extends App {
  new Control(new Version3App(), "Phase transformation")  
}

class Version3App extends Simulation {
  val defgradPlot = new Plot("Deformation Gradient")
  
  var sim: Version3Sim = _
  
  def load(c: Control) {
    c.frameTogether("Plots", defgradPlot)
    
    params.add("L", 12.0)
    params.add("dt", 0.001)
    params.add("dx", 0.01)
    params.add("omega", 1.0)
    params.add("v", 1.0)
    params.add("time")
  }
  
  def animate() {
    defgradPlot.registerLines("A0", new PointSet(0, sim.dx, sim.A0), Color.BLUE)
    defgradPlot.registerLines("A1", new PointSet(0, sim.dx, sim.A), Color.RED)
    
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
    
    sim = new Version3Sim(L=Lp, dx=dx, dt=dt, omega=omega, v=v)
    
    while (true) {

     Job.animate()
     
     val lastTime = sim.time
     while (sim.time-lastTime < 0.1) {
//       sim.laxFriedrichsStep()
//       sim.laxWendroffStep()
//       sim.semiImplicitEulerStep()
//       sim.upwindStep()
       sim.implicitUpwindStep()
     }
     
     Thread.sleep(10)
    }
  }
}