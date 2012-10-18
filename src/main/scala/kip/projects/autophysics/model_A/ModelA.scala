package kip.projects.autophysics.model_A

import kip.javasim.Random
import kip.javasim.LatticeNeighbors
import scikit.jobs.Simulation
import scikit.graphics.dim2.Grid
import scikit.jobs.Control
import scikit.jobs.Job



class Ising2d(val seed: Int, val L: Int, var T: Double) {
  var time = 0.0
  val h = 0.0
  val N  = L*L
  val spin = new Array[Int](N)
  
  val random = new Random(seed)
  
  randomize()
  
  def randomize() {
    for (i <- 0 until N)
      spin(i) = if (random.nextDouble() < 0.5) 1 else -1
  }
  
  val neighbors = {
    val nl = new LatticeNeighbors(L, L, 0.1, 1.1, LatticeNeighbors.Type.PERIODIC)
    Array.tabulate(N) { i => nl.get(i) }
  }
  
  require(neighbors.map(_.size).distinct.size == 1)
  
  def neighborSum(i: Int): Int = {
    var acc = 0
    for (j <- 0 until neighbors(i).size)
      acc += spin(neighbors(i)(j))
    acc
  }
  
  def singleStep() {
    val i = random.nextInt(N)
    val dE = 2*spin(i)*neighborSum(i) + 2*h*spin(i)
    if (dE <= 0 || (T != 0 && random.nextDouble() < math.exp(-dE/T))) {
      spin(i) = -spin(i)
    }
  }
    
  def step(mcs: Double) {
    val n = math.max(mcs * N, 0).toInt
    for (i <- 0 until n)
      singleStep();
    time += n.toDouble / N;
  }
  
  def cgField(Lp: Int): Array[Double] = {
    require(L % Lp == 0)
    val cgLength = L / Lp
    val ret = new Array[Double](Lp*Lp)
    for (xp <- 0 until Lp;
         yp <- 0 until Lp) {
      var acc = 0.0
      for (x <- xp*cgLength until (xp+1)*cgLength;
           y <- yp*cgLength until (yp+1)*cgLength) {
        acc += spin(y*L + x)
      }
      ret(yp*Lp + xp) = acc / (cgLength*cgLength)
    }
    ret
  }
}


class CGModel(val Lp: Int, val neighborDist: Double = 2.0, val maxPower: Int = 5, val updateSize: Double = 0.1) {
  var corrections = 1
  
  val neighbors = {
    val nl = new LatticeNeighbors(Lp, Lp, 0.0, neighborDist, LatticeNeighbors.Type.PERIODIC)
    Array.tabulate(Lp*Lp) { i => nl.get(i) }
  }
  
  val nneighbors = neighbors(0).size
  
  // every neighbor list should be the same size
  require(neighbors.forall(_.size == nneighbors))
  
  val alpha_p  = new Array[Double](nneighbors)
  val alpha_pp = new Array[Double](nneighbors)
  val beta = new Array[Double](maxPower)

  var phi_p: Array[Double] = new Array[Double](Lp*Lp)
  var phi_pp: Array[Double] = new Array[Double](Lp*Lp)
  var prediction: Array[Double] = new Array[Double](Lp*Lp)
  
  def addData(phi: Array[Double], update: Boolean) {
    // update model based on exact result "phi"
    if (update) {
      for (i <- 0 until Lp*Lp) {
        correct(i, phi(i))
      }
      corrections += 1
    }
    
    // store "phi" as previous state 
    phi_pp = phi_p
    phi_p = phi
    
    // make a new prediction
    for (i <- 0 until Lp*Lp) {
      prediction(i) = predict(i)
    }
  }

  def predict(i: Int): Double = {
    var acc = 0.0
    for ((j, idx) <- neighbors(i).zipWithIndex) {
      acc += alpha_p(idx) * phi_p(j)
      acc += alpha_pp(idx) * phi_pp(j)
    }
    for (k <- 2 until maxPower) {
      acc += beta(k) * math.pow(phi_p(i), k) 
    }
    acc
  }
  
  def correct(i: Int, phi_ex: Double) = {
    val phi_tilde = predict(i)
    
    val dev = phi_tilde - phi_ex
    val del = updateSize / math.sqrt(corrections)
    
    for ((j, idx) <- neighbors(i).zipWithIndex) {
      alpha_p(idx)  -= del * dev * phi_p(j)
      alpha_pp(idx) -= del * dev * phi_pp(j)
    }
    for (k <- 2 until maxPower) {
      beta(k) -= del * dev * math.pow(phi_p(i), k) 
    }
  }
  
  def printModel() {
    println("beta = "+beta.mkString("<", ",", ">"))
    println("alpha_p = "+alpha_p.mkString("<", ",", ">"))
    println("alpha_pp = "+alpha_pp.mkString("<", ",", ">"))
    println()
  }
}



object Ising2DApp {
  
  def main(args: Array[String]) {
    new Control(new Ising2DApp(), "Ising Model");
  }  
}

class Ising2DApp extends Simulation {
  val grid = new Grid("Ising spins")
  val cgGrid = new Grid("CG spins")
  val predictGrid = new Grid("CG prediction")
  
  cgGrid.setScale(-1, 1)
  predictGrid.setScale(-1, 1)
  println("set scale")
  
  var sim: Ising2d = _
  var cgModel: CGModel = _
  
  var dt: Double = _
  
  override def load(c: Control) {
    c.frame(grid, cgGrid, predictGrid)

    params.add("Seed", 69)
    params.add("L", 100)
    params.add("CG length", 4)
    params.addm("T", 2.0)
    params.addm("dt", 10.0)
  }

  override def animate() {
    val L = sim.L
    val Lp = cgModel.Lp

    sim.T = params.fget("T");
    dt = params.fget("dt");
    grid.registerData(sim.L, sim.L, sim.spin)
    
    cgGrid.registerData(Lp, Lp, sim.cgField(Lp))
    
    predictGrid.registerData(Lp, Lp, cgModel.prediction)
    
    cgModel.printModel()
  }
  
  override def clear() {
    grid.clear()
    cgGrid.clear()
    predictGrid.clear()
  }
  
  def run() {
    dt = params.fget("dt");
    val seed = params.iget("Seed");
    val L = params.iget("L");
    sim = new Ising2d(seed, L, params.fget("T"))
    
    val cg_length = params.iget("CG length")
    require(L % cg_length == 0)
    
    val Lp = L / cg_length
    
    cgModel = new CGModel(Lp, updateSize = 0.001)
    
    
    while (true) {
      
      sim.randomize()
      
      sim.step(dt)
      cgModel.addData(sim.cgField(Lp), false)
      sim.step(dt)
      cgModel.addData(sim.cgField(Lp), false)
      
      for (i <- 0 until 20) {
        sim.step(dt)
        cgModel.addData(sim.cgField(Lp), true)
        Job.animate()
      }
    }
  }
}
