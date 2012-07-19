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
  val energyPlots = new Plot("Energies")
  val stressPlot = new Plot("Stress")

  var potential = new Accumulator()
  var kinetic = new Accumulator()
  var totalEnergy = new Accumulator()
  var stressAcc = new Accumulator()
  
  var sim: SubMD2d = _
  
  def load(c: Control) {
    c.frameTogether("", canvas, energyPlots, stressPlot)
    
    params.add("ncols", 4)
    params.add("nrows", 4)
    params.add("a", 1.12246)
    params.add("dt", 0.01)
    params.addm("defgrad", 1.0)
    params.addm("target kinetic energy", 0)
  }
  
  def animate() {
    val circles = sim.p.map(p => Geom2D.circle(p.x, p.y, sim.sigma/2, Color.BLACK))
    import scala.collection.JavaConversions._
    canvas.setDrawables(circles.toSeq)
    val bds = new Bounds(0, sim.w0*sim.defgrad, 0, sim.h0, 0, 0)
    canvas.addDrawable(Geom2D.rectangle(bds, Color.GREEN))
    
    energyPlots.registerLines("Potential", potential, Color.RED)
    energyPlots.registerLines("Kinetic", kinetic, Color.BLUE)
    energyPlots.registerLines("Total", totalEnergy, Color.BLACK)
    
    stressPlot.registerLines("Stress", stressAcc, Color.BLACK)
    
    val defgrad = params.fget("defgrad")
    if (defgrad != sim.defgrad) {
      sim.applyDeformationGradient(defgrad)
      sim.initializeCrystalPositions()
      sim.applyKineticEnergy(params.fget("target kinetic energy"))
    }
  }
  
  def clear() {
    potential = new Accumulator()
    kinetic = new Accumulator()
    totalEnergy = new Accumulator()
    stressAcc = new Accumulator()
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
       
       stressAcc.accum(time, sim.convertCauchyToFirstPiolaKirchoff(sim.virialStress()))
     }
     Thread.sleep(10)
    }
  }
}
