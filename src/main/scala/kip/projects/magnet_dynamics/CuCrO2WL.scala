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

class EnergyBins(val e_lo: Double, val e_hi: Double, val n: Int) {
  val de = (e_hi-e_lo) / (n-1)
  
  def index(e: Double) = {
    require(e >= e_lo && e <= e_hi)
    math.round((e - e_lo) / de).toInt
  }
  
  def center(i: Int) = e_lo + i*de
  
  def interpolate(a: Array[Double], e: Double): Double = {
    require(e >= e_lo && e <= e_hi)
    val x = (e - e_lo) / de
    val alpha = x % 1.0
    val ilo = math.floor(x).toInt
    val ihi = math.ceil(x).toInt
    (1.0 - alpha) * a(ilo) + alpha * a(ihi) 
  }
}

trait Model {
  def energy: Double
  def newTrialChange(): Double
  def acceptTrialChange()
}

class IsingModel(L: Int, impurityDensity: Double) extends Model {
  val N = L*L
  val spin = new Array[Int](N)
  var energy: Double = _
  
  init()
  
  def init() {
    for (i <- 0 until N) {
      spin(i) = if (Math.random() < 0.5) 1 else -1
      if (Math.random() < impurityDensity)
        spin(i) = 0
    }

    energy = 0;
    for (i <- 0 until N) {
      energy += - spin(i) * sumNeighbors(i)
    }
    energy /= 2        // we double counted all interacting pairs
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
  
  var trial_i: Int = _
  var trial_de: Double = _

  def newTrialChange(): Double = {
    trial_i = (N * Math.random()).toInt
    trial_de = 2*spin(trial_i)*sumNeighbors(trial_i)
    trial_de
  }
  
  def acceptTrialChange() {
    spin(trial_i) = -spin(trial_i)
    energy += trial_de;
    trial_i = -1
    trial_de = 0.0
  }
}

class CuCrO2WL extends Simulation {
  val histogramPlot = new Plot("Histogram", "Energy", "Histogram");
  val densityPlot = new Plot("Density of states", "Energy", "Density of states");
  val heatPlot = new Plot("Heat capacity", "Temperature", "Heat capacity");

  var model: IsingModel = _
  var bins: EnergyBins = _
  var mcs: Int = _         // total MC system sweeps
  var f: Double = _        // multiplicative modification factor to g
  var iterations: Int = _  // number of reductions to f
  var g: Array[Double] = _ // logarithm of the density of states (energy argument translated by 2N)
  var H: Array[Int] = _    // histogram (reduce f when it is "flat")
  
  def load(c: Control) {
    c.frame(histogramPlot, densityPlot, heatPlot)
    params.add("L", 16)
    params.add("Impurity density", 0.0)
    params.add("Energy bins", 100)
    params.add("mcs")
    params.add("iteration")
  }

  def run() {
    model = new IsingModel(L=params.iget("L"),
                           impurityDensity=params.fget("Impurity density"))
    bins = new EnergyBins(e_lo = -2*model.N, e_hi = 2*model.N, n=(4*model.N+1)/10)
    mcs = 0
    f = Math.exp(1)
    iterations = 0
    g = Array.fill(bins.n)(0.0)
    H = Array.fill(bins.n)(0)
    
    while (true) {
      doStep();
      Job.animate();
    }
  }

  def flipSpins() {
    for (steps <- 0 until model.N) {
      val de = model.newTrialChange()
      require (model.energy+de >= bins.e_lo && model.energy+de <= bins.e_hi)
      val dg = bins.interpolate(g, model.energy) - bins.interpolate(g, model.energy+de)
      if (Math.random() < Math.exp(dg)) {
        model.acceptTrialChange()
      }
      
      val i = bins.index(model.energy)
      g(i) += Math.log(f)
      H(i) += 1
      
      Job.`yield`()
    }
  }

  def isFlat(): Boolean = {
    var netH = 0
    var numEnergies = 0
    
    for (i <- 0 until bins.n) {
      if (H(i) > 0) {
        netH += H(i)
        numEnergies += 1
      }
    }

    for (i <- 0 until bins.n)
      if (0 < H(i) && H(i) < 0.8*netH/numEnergies)
        return false

    true
  }


  def doStep() {
    val mcsMax = mcs + Math.max(200000/model.N, 1)
    while (mcs < mcsMax) {
      flipSpins()
      mcs += 1
    }
    
    if (isFlat()) {
      f = math.sqrt(f)
      iterations += 1
      H.transform(_ => 0)
    }
  }
  
  
  def animate() {
    params.set("mcs", mcs)
    params.set("iteration", iterations)

    val densityData = new DynamicArray()
    val histogramData = new DynamicArray()
    val heatData = new Accumulator(0.02)

    for (i <- 0 until bins.n) {
      if (g(i) > 0) {
        densityData.append2  (bins.center(i), g(i) - g(0))
        histogramData.append2(bins.center(i), H(i))
      }
    }
    
    for (T <- 0.5 to 5 by 0.1)
      heatData.accum(T, Thermodynamics.heatCapacity(bins, g, 1/T))
    for (T <- 1.2 to 1.8 by 0.02)
      heatData.accum(T, Thermodynamics.heatCapacity(bins, g, 1/T))
    
    
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
  def logZ(bins: EnergyBins, g: Array[Double], beta: Double): Double = {
    // m = max {e^(g - beta E)}
    var m = 0.0;
    for (i <- 0 until bins.n)
      m = Math.max(m, g(i) - beta*bins.center(i));

    //     s     = Sum {e^(g - beta E)} * e^(-m)
    // =>  s     = Z * e^(-m)
    // =>  log s = log Z - m
    var s = 0.0
    for (i <- 0 until bins.n)
      s += Math.exp(g(i) - beta*bins.center(i) - m)
    Math.log(s) + m
  }


  def heatCapacity(bins: EnergyBins, g: Array[Double], beta: Double): Double = {
    val logZ_ = logZ(bins, g, beta)

    var E_avg = 0.0;
    var E2_avg = 0.0;

    for (i <- 0 until bins.n) {
      if (g(i) != 0) {
        val e = bins.center(i)
        E_avg  += e   * Math.exp(g(i) - beta*e - logZ_)
        E2_avg += e*e * Math.exp(g(i) - beta*e - logZ_)
      }
    }

    (E2_avg - E_avg*E_avg) * beta*beta;
  }
}
