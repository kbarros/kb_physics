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
  val waveSpeed = 0.5
  def apply(strain: Double, velocity: Double, energy: Double) = waveSpeed*waveSpeed * strain
}


class SimpleShockSolver(L: Int, dx: Double, dt: Double) {
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
    val del_i = ((time * SimpleStress.waveSpeed) / dx).toInt
    Array.tabulate(L) { i =>
      val i1 = mod(i+del_i, L)
      val i2 = mod(i-del_i, L)
      0.5*(strain0(i1) + strain0(i2))
    } 
  }  
}


object SimpleShockApp extends App {
  new Control(new SimpleShockApp(), "Shock wave")  
}

class SimpleShockApp extends Simulation {
  val canvas = new Scene2D("Particles");
  val strainPlot = new Plot("Strain");

  var sim: SimpleShockSolver = _
  
  def load(c: Control) {
    c.frame(strainPlot)
    c.frame(canvas);

    params.add("L", 400)
    params.add("dx", 1.0)
    params.add("dt", 0.5)
  }
  
  def animate() {
    strainPlot.registerLines("strain1", new PointSet(0, 1, sim.sve1.strain), Color.GRAY)
    strainPlot.registerLines("strain2", new PointSet(0, 1, sim.sve2.strain), Color.BLACK)
    strainPlot.registerLines("strain exact", new PointSet(0, 1, sim.linearStrainSolution()), Color.RED)
    
    import scala.collection.JavaConversions._
    canvas.setDrawables(Seq(Geom2D.circle(1, 1, 2, Color.BLACK)))
    
//    println("Total energy " + sim.sve2.energy.sum)
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
