package kip.projects.spinglass

import scikit.graphics.dim2.Grid
import scikit.jobs.Control
import scikit.jobs.Job
import scikit.jobs.Simulation
import scikit.jobs.params.DoubleValue
import scikit.jobs.params.Parameters
import scikit.numerics.Math2;
import scala.util.Random


class Annni2dSim(params: Parameters) {
  val random = new Random();

  val L = params.iget("L");
  val N = L*L;
  val occupiedFraction = params.fget("Occupied fraction");
  
  val spin = new Array[Int](N)
  val occupied = new Array[Int](N) // 1 if occupied, 0 if vacant
  val field = new Array[Double](N)
  
  random.setSeed(params.iget("Random seed"));
  for (i <- 0 until N) {
    spin(i) = if (random.nextDouble() < 0.5) 1 else -1
    field(i) = random.nextGaussian() 
    occupied(i) = if (random.nextDouble() < occupiedFraction) 1 else 0
  }

  var T: Double = _
  var fieldMean: Double = _
  var fieldStrength: Double = _
  readParams()
    
  def readParams() {
    T = params.fget("T")
    fieldMean = params.fget("Field mean")
    fieldStrength = params.fget("Field std. dev.")
  }
    
  def buildNeighbors(i: Int): (Array[Int], Array[Double]) = {
    val displacements: Seq[(Int, Int, Double)] = Seq(
      // ferromagnetic nn
      (+1, 0, -1),
      (0, +1, -1),
      (-1, 0, -1),
      (0, -1, -1),
      
      // anti-ferrogment nnn
      (+1, +1, +1),
      (-1, +1, +1),
      (+1, -1, +1),
      (-1, -1, +1)
    )
    
    val y = i/L
    val x = i%L
    val (idx, strength) = (for ((dx, dy, a) <- displacements) yield {
      val xp = (x + dx + L) % L
      val yp = (y + dy + L) % L;
      (yp * L + xp, a)
    }).unzip
    (idx.toArray, strength.toArray)
  }
    
  val neighbors = Array.tabulate[(Array[Int], Array[Double])](N)(buildNeighbors _)
    
  // energy contribution between spin i and neighbors
  def neighborEnergy(i: Int): Double = { 
    val (n, strength) = neighbors(i)
    var acc = 0d
    for (nidx <- 0 until n.size) {
      val j = n(nidx)
      val a = strength(nidx)
      acc += occupied(i)*occupied(j)*a*spin(i)*spin(j)
    }
    acc
  }
    
  // energy contribution between spin i and external field
  def fieldEnergy(i: Int): Double = {
    - (occupied(i) * spin(i)) * (fieldMean + fieldStrength * field(i))
  }
    
  def netEnergy(): Double = {
    var acc = 0d
    for (i <- 0 until N) {
      acc += neighborEnergy(i) / 2.0 // double counting
      acc += fieldEnergy(i)
    }
    acc
  }
  
  def step() {
    val i = random.nextInt(N);
    
    val e0 = neighborEnergy(i) + fieldEnergy(i)
    spin(i) *= -1
    val e1 = neighborEnergy(i) + fieldEnergy(i)
    spin(i) *= -1
    
    val de = e1 - e0
    if (de <= 0 || (T != 0 && random.nextDouble() < math.exp(-de/T))) {
      spin(i) *= -1
    }
  }
}


object Annni2d extends App {
  new Control(new Annni2d(), "Strain glass");  
}

class Annni2d extends Simulation {
  val grid = new Grid("Ising spins");
  var sim: Annni2dSim = _

  def load(c: Control) {
    grid.setScale(-1, +1);
    c.frame(grid);
    params.add("L", 32);
    params.add("Random seed", 0);
    params.add("Occupied fraction", 1.0);
    
    params.addm("T",               new DoubleValue(0, 0, 10).withSlider());
    params.addm("Field mean",      new DoubleValue(0, -10, 10).withSlider())
    params.addm("Field std. dev.", new DoubleValue(0, 0, 2).withSlider())
  }

  def animate() {
    sim.readParams();
    grid.registerData(sim.L, sim.L, sim.spin);
  }

  def clear() {
    grid.clear();
  }

  def run() {
    sim = new Annni2dSim(params);    

    while (true) {
      sim.step();           
      Job.animate();
    }
  }
}
