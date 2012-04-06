package kip.projects.langevin

import scikit.util.Utilities.format

import java.awt.Color

import scikit.dataset.Histogram
import scikit.graphics.dim2.Grid
import scikit.graphics.dim2.Plot
import scikit.jobs.Control
import scikit.jobs.Job
import scikit.jobs.Simulation
import scikit.jobs.params.ChoiceValue

object TDGLCollectorApp extends App {
  new Control(new TDGLCollectorApp(), "Allen Cahn Coarsening")
}

class TDGLCollectorApp extends Simulation {
  val timesPlot = new Plot("Times")
  var times = new Histogram(1.0)
  var sim: TDGL = _
  val dimensions = 3

  def load(c: Control) {
    c.frame(timesPlot)
    params.addm("dt", 5.0)
    params.add("L", 25.0)
    params.add("Random seed", 0)
    params.add("dx", 1.0)
    params.add("r")
    params.add("Initial cond", new ChoiceValue("Small", "Small w/ shift", "Ising"))
    params.add("Time/L^2")
    params.add("Energy")
    params.add("Stripe fraction")

    params.set("r", 1.0)
  }

  def animate() {
    timesPlot.registerBars("Escape times", times, Color.BLUE)

    params.set("Time/L^2", format(sim.t))
    params.set("Energy", format(sim.freeEnergy))
  }

  def clear() {
  }

  def run() {
    times = new Histogram(0.1)
    sim = new TDGL(params, dimensions)

    val initialType = params.sget("Initial cond")

    var slabCount = 0
    for (iter <- 1 until Int.MaxValue) {
      params.set("Random seed", iter)
      sim.random.setSeed(iter)
      
      if ("Small".equals(initialType))
        sim.randomize()
      else if ("Small w/ shift".equals(initialType))
        sim.randomizeAndShift()
      else if ("Ising".equals(initialType))
        sim.randomizeIsing()
      else
        System.out.println("Unknown initial condition")

      Job.animate()

      val b = coarsensToStripe1
      if (b) {
        slabCount += 1
        val stripeEnergy = (2 * TDGL.surfaceEnergyDensity * Seq.fill(dimensions-1)(sim.L).product)
        val ratio = sim.freeEnergy / stripeEnergy
        if (ratio > 0.667)
          println(iter + " " + ratio)
      }

      val fraction = slabCount.toDouble / iter
      val plusMinus = math.sqrt(slabCount) / iter
      params.set("Stripe fraction", slabCount+"/"+iter+" = "+format(fraction)+" +- "+format(plusMinus))
      Job.animate()
    }
  }

  def coarsensToStripe1: Boolean = {
    val L2 = sim.L * sim.L
    val maxTime = 8.0 * L2 // L ~ t^{1/2}
    val stripeEnergy = (2 * TDGL.surfaceEnergyDensity * Seq.fill(dimensions-1)(sim.L).product)

    while (sim.t < maxTime && sim.freeEnergy > stripeEnergy / 2) {
      for (i <- 0 until 5)
        sim.simulate()
      Job.animate()
    }

    times.accum(sim.t / L2)

    sim.freeEnergy > stripeEnergy / 2
  }
}
