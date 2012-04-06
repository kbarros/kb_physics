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
  
  val spin = new Array[Int](N)
  val occupied = new Array[Int](N) // occupation priority indexed by spin
  val field = new Array[Double](N)
  
  random.setSeed(params.iget("Random seed"));
  for (i <- 0 until N) {
    spin(i) = if (random.nextDouble() < 0.5) 1 else -1
    field(i) = random.nextGaussian() 
    occupied(i) = i
  }
  for (i <- 0 until N-1) {
    val k = i + random.nextInt(N - i)
    // swap occupied(i) and occupied(k)
    val t = occupied(i)
    occupied(i) = occupied(k)
    occupied(k) = t
  }
  
  var occupiedCnt: Int = _
  var T: Double = _
  var fieldMean: Double = _
  var fieldStrength: Double = _
  var J = Array[Double](-1, 1)
  
  readParams()
  
  def readParams() {
    T = params.fget("T")
    fieldMean = params.fget("Field mean")
    fieldStrength = params.fget("Field std. dev.")
    J(1) = params.fget("J2 / J1")
    occupiedCnt = (params.fget("Occupied fraction") * N).toInt max 0 min (N-1)
    
    for (i <- 0 until N) {
      if (occupied(i) <= occupiedCnt && spin(i) == 0)
        spin(i) = if (random.nextDouble() < 0.5) 1 else -1
      else if (occupied(i) > occupiedCnt && spin(i) != 0)
        spin(i) = 0
    }
  }
  
  def buildNeighbors(i: Int): (Array[Int], Array[Int]) = {
    val displacements: Seq[(Int, Int, Int)] = Seq(
      // ferromagnetic nn
      (+1, 0, 0),
      (0, +1, 0),
      (-1, 0, 0),
      (0, -1, 0),
      
//      // anti-ferrogment nnn
//      (+1, +1, 1),
//      (-1, +1, 1),
//      (+1, -1, 1),
//      (-1, -1, 1)
      
      // anti-ferrogment nnn
      (+2, 0, 1),
      (-2, 0, 1)
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
  
  val neighbors = Array.tabulate[(Array[Int], Array[Int])](N)(buildNeighbors _)
    
  // energy contribution between spin i and neighbors
  def neighborEnergy(i: Int): Double = { 
    val (n, ji) = neighbors(i)
    var acc = 0d
    for (nidx <- 0 until n.size) {
      val j = n(nidx)
      val a = J(ji(nidx))
      acc += a*spin(i)*spin(j)
    }
    acc
  }
    
  // energy contribution between spin i and external field
  def fieldEnergy(i: Int): Double = {
    - spin(i) * (fieldMean + fieldStrength * field(i))
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
    params.add("L", 64);
    params.add("Random seed", 0);
    
    params.addm("T",               new DoubleValue(1, 0, 2).withSlider());
    params.addm("J2 / J1",         new DoubleValue(0.5, 0.1, 1).withSlider())
    params.addm("Field mean",      new DoubleValue(0, -10, 10).withSlider())
    params.addm("Field std. dev.", new DoubleValue(0, 0, 2).withSlider())
    params.addm("Occupied fraction", new DoubleValue(1, 0, 1).withSlider())
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
      for (i <- 0 until 1000)
        sim.step();
      Job.animate();
    }
  }
}
