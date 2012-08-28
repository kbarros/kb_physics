package kip.projects.spinglass.v2

import scala.util.Random
import kip.math.fft.FFTReal
import scikit.jobs.Control
import scikit.jobs.Simulation
import scikit.graphics.dim2.Grid
import scikit.jobs.Job

class LongRangeDisorderSim(seed: Int, val L: Int, dx: Double, var dt: Double, alpha: Double, J: Double,
                           var H: Double, var dilution: Double, var T: Double) {
  val rand = new Random(seed)

  val N = L*L
  
  val phi = new Array[Double](N)
  val phibar = new Array[Double](N) 
  val efield = new Array[Double](N) 
  phi(0) = 1

  val interaction = new Array[Double](N)
  
  buildInteraction()
  println("inter min " + interaction.min + " max " + interaction.max)
  
  val H_disorder = Array.tabulate[Double](N) { _ => rand.nextGaussian() }
  val site_disorder = Array.tabulate[Double](N) { _ => rand.nextDouble() }
  
  val fft = new FFTReal(dim=Array(L, L))
  
  def buildInteraction() {
    // long range repulsion
    for (y <- 0 until L; x <- 0 until L) {
      val i = x + L*y
      val (r, theta) = {
        val xp = if (x < L/2) x else x - L
        val yp = if (y < L/2) y else y - L
        val r = dx * math.sqrt((xp*xp + yp*yp)) 
        val theta = math.atan2(yp, xp)
        (r, theta)
      }
      
      val ferroRange = 1.1 * dx
      interaction(i) = if (r < ferroRange) -J else math.pow(r, -alpha)  * math.cos(4*theta)
    }
  }
  
  // \sum phi phibar / 2 + phi^4 / 4 + H phi
  def energy(): Double = {
    fft.convolve(phi, interaction, phibar)
    for (i <- 0 until N) {
      efield(i) = (phi(i)*phibar(i)/2 + phi(i)*phi(i)*phi(i)*phi(i)/4 + H*H_disorder(i)*phi(i))
    }
    efield.sum / N
  }

  def step() {
    fft.convolve(phi, interaction, phibar)
    
    for (i <- 0 until N) {
      phi(i) -= dt * (phibar(i) + phi(i)*phi(i)*phi(i) + H*H_disorder(i)) + math.sqrt(dt*2*T)*rand.nextGaussian()
      
      if (site_disorder(i) < dilution)
        phi(i) = 0
    }
  }
}

object LongRangeDisorderApp extends App {
  new Control(new LongRangeDisorderApp(), "Long Range Disorder App");  
}

class LongRangeDisorderApp extends Simulation {
  val grid1 = new Grid("Phi")
  val grid2 = new Grid("Energy")
  var sim: LongRangeDisorderSim = _
  
  def load(c: Control) {
    c.frame(grid1, grid2)
    params.add("seed", 0)
    params.add("L", 50)
    params.add("dx", 1.0)
    params.add("alpha", 2.0)
    params.add("J", 2.0)
    params.addm("dt", 0.005)
    params.addm("H", 0.0)
    params.addm("dilution", 0.0)
    params.addm("T", 10.0)
    
    params.add("energy")
  }

  def animate() {
    grid1.registerData(sim.L, sim.L, sim.phi)
    grid2.registerData(sim.L, sim.L, sim.efield)
    
    params.set("energy", sim.energy())
    
    sim.dt = params.fget("dt")
    sim.H = params.fget("H")
    sim.dilution = params.fget("dilution")
    sim.T = params.fget("T")
  }
  
  def clear() {
    grid1.clear()
    grid2.clear()
  }

  def run() {
    val seed = params.iget("seed")
    val L = params.iget("L")
    val dx = params.fget("dx")
    val dt = params.fget("dt")
    val alpha = params.fget("alpha")
    val J = params.fget("J")
    val H = params.fget("H") 
    val dilution = params.fget("dilution") 
    val T = params.fget("T")
    sim = new LongRangeDisorderSim(seed=seed, L=L, dx=dx, dt=dt, alpha=alpha,
                                   J=J, H=H, dilution=dilution, T=T)
   
    while (true) {
      Job.animate();
      for (i <- 0 until 1)
        sim.step()
    }
  }
}
