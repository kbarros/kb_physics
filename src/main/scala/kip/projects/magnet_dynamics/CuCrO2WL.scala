package kip.projects.magnet_dynamics


import java.awt.Color;

import scikit.dataset.Accumulator;
import scikit.dataset.DynamicArray;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;

object CuCrO2WL {
  def main (args: Array[String]) {
    new Control(new CuCrO2WL(), "Wang Landau of CuCrO2")
  }
}

trait Model {
  def energy()
  def newTrialChange(): Double
  def applyChange()
}

class IsingModel(L: Int) extends Model {
  def applyChange(): Unit = ???
  def energy(): Unit = ???
  def newTrialChange(): Double = ???
}

class CuCrO2WL extends Simulation {
  val histogramPlot = new Plot("Histogram", "Energy", "Histogram");
  val densityPlot = new Plot("Density of states", "Energy", "Density of states");
  val heatPlot = new Plot("Heat capacity", "Temperature", "Heat capacity");

  var mcs: Int = _
  var L: Int = _
  var N: Int = _
  var density: Double = _  // percentage of (spin 0) magnetic impurities
  var g: Array[Double] = _ // logarithm of the density of states (energy argument translated by 2N)
  var H: Array[Int] = _    // histogram (reduce f when it is "flat")
  var E: Int = _           // energy of current spin configuration (translated by 2N)
  var spin: Array[Int] = _
  var f: Double = _        // multiplicative modification factor to g
  var iterations: Int = _  // number of reductions to f
  
  def load(c: Control) {
    c.frame(histogramPlot, densityPlot, heatPlot)
    params.add("L", 16)
    params.add("Impurity density", 0.2)
    params.add("mcs")
    params.add("iteration")
  }

  def run() {
    L = params.iget("L")
    density = params.fget("Impurity density")

    mcs = 0
    N = L*L
    f = Math.exp(1)
    iterations = 0

    spin = new Array[Int](N)
    for (i <- 0 until N) {
      spin(i) = if (Math.random() < 0.5) 1 else -1
      if (Math.random() < density)
        spin(i) = 0
    }

    g = new Array[Double](4*N + 1)
    H = new Array[Int]   (4*N + 1)
    for (e <- 0 to 4*N) {
      g(e) = 0
      H(e) = 0
    }
    
    E = 0;
    for (i <- 0 until N) {
      E += - spin(i) * sumNeighbors(i)
    }
    E /= 2        // we double counted all interacting pairs
    E += 2*N      // translate energy by 2*N to facilitate array access

    while (true) {
      doStep();
      Job.animate();
    }
  }


  def sumNeighbors(i: Int): Int = {
    var u = i - L;
    var d = i + L;
    var l = i - 1;
    var r = i + 1;

    if (u < 0)        u += N;
    if (d >= N)       d -= N;
    if (i % L == 0)   l += L;
    if (r % L == 0)   r -= L;
    spin(u) + spin(d) + spin(l) + spin(r)
  }

  def flipSpins() {
    for (steps <- 0 until N) {
      val i = (N * Math.random()).toInt
      val dE = 2*spin(i)*sumNeighbors(i)

      if (Math.random() < Math.exp(g(E) - g(E + dE))) {
        spin(i) = -spin(i)
        E += dE;
      }

      g(E) += Math.log(f)
      H(E) += 1
      
      Job.`yield`()
    }
  }

  def isFlat(): Boolean = {
    var netH = 0
    var numEnergies = 0
    
    for (e <- 0 to 4*N) {
      if (H(e) > 0) {
        netH += H(e)
        numEnergies += 1
      }
    }

    for (e <- 0 to 4*N)
      if (0 < H(e) && H(e) < 0.8*netH/numEnergies)
        return false

    true
  }


  def doStep() {
    val mcsMax = mcs + Math.max(100000/N, 1)
    while (mcs < mcsMax) {
      flipSpins()
      mcs += 1
    }
    
    if (isFlat()) {
      f = Math.sqrt(f)
      iterations += 1
      for (e <- 0 to 4*N)
        H(e) = 0
    }
  }
  
  
  def animate() {
    params.set("mcs", mcs)
    params.set("iteration", iterations)

    val densityData = new DynamicArray()
    val histogramData = new DynamicArray()
    val heatData = new Accumulator(0.02)

    for (e <- 0 to 4*N) {
      if (g(e) > 0) {
        densityData.append2  (e - 2*N, g(e) - g(0))
        histogramData.append2(e - 2*N, H(e))
      }
    }
    
    for (T <- 0.5 to 5 by 0.1)
      heatData.accum(T, Thermodynamics.heatCapacity(N, g, 1/T))
    for (T <- 1.2 to 1.8 by 0.02)
      heatData.accum(T, Thermodynamics.heatCapacity(N, g, 1/T))
    
    
    densityPlot.setAutoScale(true)
    densityPlot.registerPoints("Density", densityData, Color.BLUE)
    
    histogramPlot.registerPoints("Histogram", histogramData, Color.BLACK)
    
    heatPlot.setAutoScale(true)
    heatPlot.registerLines("Heat capacity", heatData, Color.RED)
  }
  
  def clear() {
    densityPlot.clear();
    histogramPlot.clear();
    heatPlot.clear();
  }
}

object Thermodynamics {
  def logZ(N: Int, g: Array[Double], beta: Double): Double = {
    // m = max {e^(g - beta E)}
    var m = 0.0;
    for (E <- -2*N to 2*N)
      m = Math.max(m, g(E+2*N) - beta*E);

    //     s     = Sum {e^(g - beta E)} * e^(-m)
    // =>  s     = Z * e^(-m)
    // =>  log s = log Z - m
    var s = 0.0
    for (E <- -2*N to 2*N)
      s += Math.exp(g(E+2*N) - beta*E - m)
    Math.log(s) + m
  }


  def heatCapacity(N: Int, g: Array[Double], beta: Double): Double = {
    val logZ_ = logZ(N, g, beta)

    var E_avg = 0.0;
    var E2_avg = 0.0;

    for (E <- -2*N to 2*N) {
      if (g(E+2*N) != 0) {
        E_avg  += E   * Math.exp(g(E+2*N) - beta*E - logZ_)
        E2_avg += E*E * Math.exp(g(E+2*N) - beta*E - logZ_)
      }
    }

    (E2_avg - E_avg*E_avg) * beta*beta;
  }
}
