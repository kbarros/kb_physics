package kip.projects.dib_nano

import java.awt.Color
import scala.math._
import scikit.dataset.PointSet
import scikit.graphics.dim2.Plot
import scikit.jobs.Control
import scikit.jobs.Job
import scikit.jobs.Simulation
import scikit.dataset.Accumulator
import kip.javasim.Random



object DibNanoApp extends App {
  new Control(new DibNanoApp(), "Shock wave")  
}

class DibNanoApp extends Simulation {
  val phiPlot = new Plot("Phi")
  val phitPlot = new Plot("Phi_t")
  
  val phiAcc = new Accumulator()
  val phitAcc = new Accumulator()
  
  val rand = new Random()
  
  def load(c: Control) {
    c.frameTogether("Plots", phiPlot, phitPlot)

    params.addm("dt", 0.01)
    params.addm("i_B", 1.0)
    params.addm("T", 0.1)
    params.addm("r_s", 1.0)
    params.addm("tau", 1.0)
    params.addm("beta", 1.0)
    params.add("phi0", 0.0)
    params.add("phit0", 0.0)
    params.add("eta0", 0.0)
  }
  
  def animate() {
    phiPlot.registerLines("phi", phiAcc, Color.RED)
    phitPlot.registerLines("phi_t", phitAcc, Color.RED)
  }
  
  def clear() {
    phiPlot.clear()
    phitPlot.clear()
    
    phiAcc.clear()
    phitAcc.clear()
  }
  
  def run() {
    val phi0 = params.fget("phi0")
    val phit0 = params.fget("phit0")
    val eta0 = params.fget("eta0")
    
    var t = 0.0
    var a = 0.0
    var phi = phi0
    var phit = phit0
    var eta = eta0
    
    while (true) {
      val dt = params.fget("dt")
      val i_B = params.fget("i_B")
      val T = params.fget("T")
      val r_s = params.fget("r_s")
      val tau = params.fget("tau")
      val beta = params.fget("beta")
      
      for (i <- 0 until 100) {
        val D = 2*T*beta/r_s
        val f = -beta*phit - (beta/(tau*r_s)) * a - abs(sin(phi)) + i_B
        phit += dt*f + dt*eta + sqrt(2*T*beta*dt) * rand.nextGaussian()
        phi  += dt*phit
        eta  += - (dt/tau) * eta + (sqrt(2*D*dt)/tau) * rand.nextGaussian()
        a    += - (dt/tau) * a + dt * phit
        t    += dt
      }
      
      phiAcc.accum(t, phi)
      phitAcc.accum(t, phit)
      
      Job.animate()

      Thread.sleep(10)
    }
  }
}
