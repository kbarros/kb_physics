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
import scikit.graphics.dim2.Plot




// Stress as a function of strain and energy
trait StressFn {
  def apply(energy: Double, strain: Double): Double
}

object SimpleStress extends StressFn {
  val speed = 0.5
  def apply(energy: Double, strain: Double) = speed*speed * strain
}



// A system of equations
//   d_t w + d_x f = 0
// where w and f(w) are vectors representing conserved fields and their fluxes
//
trait ConservationLaws {
  val L: Int
  val ws: Seq[Array[Double]]
  def fs(ws: Seq[Array[Double]]): Seq[Array[Double]]
}


trait CLIntegrator {  
  def halfIter(cl: ConservationLaws, dx: Double, dt: Double)
  
  // Integrate conserved fields cl.ws by timestep dt
  def iter(cl: ConservationLaws, dx: Double, dt: Double) {
    // two half steps on a staggered grid
    halfIter(cl, dx, 0.5*dt)
    halfIter(cl, dx, 0.5*dt)
    
    // rotate back to original index <-> coordinate map
    for (w <- cl.ws) {
      val wn = w(w.size-1)
      for (i <- w.size-1 to 1 by -1) w(i) = w(i-1)
      w(0) = wn
    }
  }
}

object CLIntegratorFirstOrder extends CLIntegrator {  
  def halfIter(cl: ConservationLaws, dx: Double, dt: Double) {
    val fs = cl.fs(cl.ws)
    for ((w, f) <- cl.ws zip fs) {
      val wc = w.clone()
      for (i <- 0 until cl.L) {
        val ip = (i+1) % cl.L
        w(i) = 0.5 * (wc(i) + wc(ip)) - (dt / dx) * (f(ip) - f(i))
      }
    }
  }
}

object CLIntegratorSecondOrder extends CLIntegrator {  
  def minmod(x: Double, y: Double): Double = {
    if (x > 0 == y > 0) 0.5*(x+y)
    else 0
  }
  
  def derivative(w: Array[Double], dx: Double): Array[Double] = {
    val L = w.size
    Array.tabulate(L) { i =>
      val im = (i-1+L)%L
      val ip = (i+1)%L
      minmod(w(ip)-w(i), w(i)-w(im)) / dx 
    }
  }
  
  def halfIter(cl: ConservationLaws, dx: Double, dt: Double) {
    val L = cl.L
    val fs = cl.fs(cl.ws)
    
    // calculate w^(n+1/2)
    val wsp = for ((w, f) <- cl.ws zip fs) yield {
      val df = derivative(f, dx)
      Array.tabulate(L)(i => w(i) - (dt/(2*dx)) * df(i))
    }
    
    // calculate f^(n+1/2)
    val fsp = cl.fs(wsp)
    
    // calculate w^(n+1)
    for (k <- 0 until cl.ws.size) {
      val w  = cl.ws(k)
      val wc = w.clone()
      val dw = derivative(wc, dx)
      val fp = fsp(k)
      
      for (i <- 0 until L) {
        val ip = (i+1)%L
        w(i) = 0.5*(wc(i) + wc(ip)) + (dx/8) * (dw(i) - dw(ip)) - (dt/dx) * (fp(ip) - fp(i))
      }
    }
  }
}


class StrainVelocityEnergy(val L: Int, strain0: Array[Double], sf: StressFn) extends ConservationLaws {
  val strain: Array[Double] = strain0.clone()
  val velocity: Array[Double] = new Array[Double](L)
  val energy: Array[Double] = new Array[Double](L)
    
  val ws = Seq(strain, velocity, energy)
  
  def fs(ws: Seq[Array[Double]]): Seq[Array[Double]] = {
    val w_strain = ws(0)
    val w_velocity = ws(1)
    val w_energy = ws(2)
    
    val stress = Array.tabulate(L)(i => sf(energy=w_energy(i), strain=w_strain(i))) 
      
    val sflux = Array.tabulate(L)(i => -w_velocity(i))
    val vflux = Array.tabulate(L)(i => -stress(i))
    val eflux = Array.tabulate(L)(i => stress(i) * w_velocity(i))
    
    Seq(sflux, vflux, eflux)
  }
  
}


class ShockSolver(L: Int, dx: Double, dt: Double) {
  // initial strain conditions
  def strain0 = Array.tabulate(L)(i => if (i < L/2) 0.01 else 0.0)
  
  var sve1 = new StrainVelocityEnergy(L, strain0, SimpleStress)  
  var sve2 = new StrainVelocityEnergy(L, strain0, SimpleStress)  
  var time = 0.0
  
  def step() {
    CLIntegratorFirstOrder.iter(sve1, dx, dt)
    CLIntegratorSecondOrder.iter(sve2, dx, dt)
    time += dt
  }
  
  def linearStrainSolution(): Array[Double] = {
    def mod(x: Int, y: Int): Int = {
      val rem = x % y
      if (rem < 0) rem + y else rem
    }
    val del_i = ((time * SimpleStress.speed) / dx).toInt
    Array.tabulate(L) { i =>
      val i1 = mod(i+del_i, L)
      val i2 = mod(i-del_i, L)
      0.5*(strain0(i1) + strain0(i2))
    } 
  }  
}


object ShockPropagation extends App {
  new Control(new ShockPropagation(), "Shock wave")  
}

class ShockPropagation extends Simulation {
  val strainPlot = new Plot("Strain");

  var sim: ShockSolver = _
  
  def load(c: Control) {
    c.frame(strainPlot)
    
    params.add("L", 400)
    params.add("dx", 1.0)
    params.add("dt", 0.5)
  }
  
  def animate() {
    strainPlot.registerLines("strain1", new PointSet(0, 1, sim.sve1.strain), Color.GRAY)
    strainPlot.registerLines("strain2", new PointSet(0, 1, sim.sve2.strain), Color.BLACK)
    strainPlot.registerLines("strain exact", new PointSet(0, 1, sim.linearStrainSolution()), Color.RED)
    
//    println("Total energy " + sim.sve2.energy.sum)
  }

  def clear() {
  }
  
  def run() {
    val L = params.iget("L")
    val dx = params.fget("dx")
    val dt = params.fget("dt")
    sim = new ShockSolver(L=L, dx=dx, dt=dt)
    
    while (true) {
     Job.animate()
     for (i <- 0 until 1)
       sim.step()
     Thread.sleep(10)
    }
  }
}
