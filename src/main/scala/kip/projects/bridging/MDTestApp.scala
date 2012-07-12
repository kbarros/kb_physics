package kip.projects.bridging

import scikit.dataset.PointSet
import scikit.jobs.Control
import scikit.jobs.Job
import scikit.jobs.Simulation
import scikit.jobs.params.DoubleValue
import scikit.jobs.params.Parameters
import scikit.numerics.Math2
import scala.util.Random
import scala.math.sqrt
import java.awt.Color
import scikit.graphics.dim2.{Scene2D, Geom2D, Plot}
import scikit.dataset.Accumulator
import scikit.util.Bounds



object MDTestApp extends App {
  new Control(new MDTestApp(), "Shock wave")  
}

class MDTestApp extends Simulation {
  val canvas = new Scene2D("Particles")
  val plots = new Plot("Energies")
  
  var potential = new Accumulator()
  var kinetic = new Accumulator()
  var totalEnergy = new Accumulator()
  
  var sim: SubMD2d = _
  
  def load(c: Control) {
    c.frame(canvas)
    c.frame(plots)
    
    params.add("ncols", 4)
    params.add("nrows", 4)
    params.add("a", 1.12246)
    params.add("dt", 0.01)
    params.addm("strain", 1.0)
    params.addm("target energy", 0)
  }
  
  def animate() {
    val circles = sim.p.map(p => Geom2D.circle(p.x, p.y, sim.sigma/2, Color.BLACK))
    import scala.collection.JavaConversions._
    canvas.setDrawables(circles.toSeq)
    val bds = new Bounds(0, sim.wu*sim.strain, 0, sim.hu, 0, 0)
    canvas.addDrawable(Geom2D.rectangle(bds, Color.GREEN))
    
    plots.registerLines("Potential", potential, Color.RED)
    plots.registerLines("Kinetic", kinetic, Color.BLUE)
    plots.registerLines("Total", totalEnergy, Color.BLACK)
    
    val strain = params.fget("strain")
    if (strain != sim.strain) {
      sim.applyStrain(strain)
      sim.initializeCrystalPositions()
      sim.applyEnergy(params.fget("target energy"))
    }
  }
  
  def clear() {
    potential = new Accumulator()
    kinetic = new Accumulator()
    totalEnergy = new Accumulator()
  }
  
  def run() {
    val ncols = params.iget("ncols")
    val nrows = params.iget("nrows")
    val a = params.fget("a")
    val dt = params.fget("dt")
    
    sim = new SubMD2d(ncols=ncols, nrows=nrows, a=a, dt=dt)
    var time = 0d
    
    while (true) {
     Job.animate()
     for (i <- 0 until 1) {
       sim.verletStep()
       time += sim.dt
       
       val ke = sim.kineticEnergy()
       val pe = sim.potentialEnergy()
       potential.accum(time, pe)
       kinetic.accum(time, ke)
       totalEnergy.accum(time, pe+ke)
     }
     Thread.sleep(10)
    }
  }
}
