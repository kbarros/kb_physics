package kip.projects.kondo_rkky

import kip.graphics.Bounds3d
import kip.graphics.RetainedScene
import kip.math.Vec3
import scikit.graphics.dim2.Grid
import scikit.jobs.Control
import scikit.jobs.Job
import scikit.jobs.Simulation
import scikit.graphics.dim3.Grid3D
import scikit.jobs.params.DoubleValue
import scala.util.Random


class HeisenbergNNSim(val L: Int) {
  val rand = new Random(0)
  
  val J1 = 1
  
  var T: Double = _
  var H: Double = _
  var J4: Double = _
  var anisotropy: Double = _
  var dS: Double = _
  
  var acceptSteps = 1
  var rejectSteps = 1
  
  val lat = new Lattice3d(L, L, L)

  val sx = new Array[Double](L*L*L)
  val sy = new Array[Double](L*L*L)
  val sz = new Array[Double](L*L*L)
  randomizeSpins()
//  ferromagnetizeSpins()
  
  def ferromagnetizeSpins() {
    for (i <- 0 until L*L*L) {
      sx(i) = 1
      sy(i) = 0
      sz(i) = 0
    }    
  }
  
  def randomizeSpins() {
    for (i <- 0 until L*L*L) {
      sx(i) = rand.nextGaussian()
      sy(i) = rand.nextGaussian()
      sz(i) = rand.nextGaussian()
      normalizeSpin(i)
    }
  }
  
  // S_i dot S_j
  def bond(i: Int, j:Int) = {
    sx(i)*sx(j) + sy(i)*sy(j) + sz(i)*sz(j)
  }
  
  def spinNorm(i: Int): Double = {
    math.sqrt(sx(i)*sx(i) + sy(i)*sy(i) + sz(i)*sz(i))
  }
  
  def normalizeSpin(i: Int) {
    val norm = spinNorm(i)
    sx(i) /= norm
    sy(i) /= norm
    sz(i) /= norm
  }
  
  def spinEnergy(i: Int): Double = {
    var acc = 0.0
    for (d <- 0 until 6) {
      val j = lat.periodicNeighbors(i)(d)
      val k = lat.periodicNeighbors(j)(d)
      acc += - J1 * bond(i, j) - J4 * bond(i, k)
    }
    
    acc += - H * sz(i)
    acc += - anisotropy * sz(i) * sz(i)
    
    acc
  }
  
  def energy(): Double = {
    var acc = 0.0
    for (i <- 0 until L*L*L)
      acc += spinEnergy(i)
    acc
  }
  
  def step() {
    for (i <- 0 until L*L*L)
      singleStep()
  }
  
  def singleStep() {
    val i = rand.nextInt(L*L*L)
    
    val e1 = spinEnergy(i)
    
    val oldSx = sx(i)
    val oldSy = sy(i)
    val oldSz = sz(i)
    
    sx(i) += dS * (rand.nextDouble() - 0.5)
    sy(i) += dS * (rand.nextDouble() - 0.5)
    sz(i) += dS * (rand.nextDouble() - 0.5)
    normalizeSpin(i)
    
    val e2 = spinEnergy(i)
    
    val de = e2 - e1
    if (de < 0 || rand.nextDouble() < math.exp(-de/T)) {
      // accept
      acceptSteps += 1
    }
    else {
      // reject
      sx(i) = oldSx
      sy(i) = oldSy
      sz(i) = oldSz
      rejectSteps += 1
    }
  }
}




object HeisenbergNNApp extends App {
  new Control(new HeisenbergNNApp(), "Heisenberg next-nearest neighbor model")  
}

class HeisenbergNNApp extends Simulation {
  val grid3d = new Grid3D("S_x")
  var sim: HeisenbergNNSim = _
  val rs = new RetainedScene(Bounds3d(Vec3(0,0,0),Vec3(1,1,0)), sizew=600, sizeh=600, cameraDistance=1.2)
  
  def load(c: Control) {
    c.frame(grid3d)
    params.add("L", 10)
    params.addm("T", 0.0)
    params.addm("J4", 0.0)
    params.addm("H", 0.0)
    params.addm("anisotropy", 0.0)
    params.addm("dS", 0.5)
    params.addm("slice", new DoubleValue(0, 0, 1).withSlider)
    params.add("energy")
    params.add("accept ratio")
    params.add("mc steps")
  }

  def setSimParams() {
    sim.T = params.fget("T")
    sim.J4 = params.fget("J4")
    sim.H = params.fget("H")
    sim.anisotropy = params.fget("anisotropy")
    sim.dS = params.fget("dS")    
  }
  
  def animate() {
    setSimParams()
    val slice = math.min((params.fget("slice") * sim.L).toInt, sim.L-1)
    
    val field = new Array[Double](3*sim.L*sim.L)
    for (i <- 0 until sim.L*sim.L) {
      val ip = i + slice*sim.L*sim.L
      field(0 + 3*i) = sim.sx(ip)
      field(1 + 3*i) = sim.sy(ip)
      field(2 + 3*i) = sim.sz(ip)
    }
    val viz = new SpinViz(sim.L, sim.L)
    viz.drawSpins(field, rs)
    
    grid3d.setScale(-1, 1)
    grid3d.registerData(sim.L, sim.L, sim.L, sim.sz.map(- _))
    
    params.set("energy", sim.energy())
    params.set("accept ratio", sim.acceptSteps.toDouble / (sim.acceptSteps + sim.rejectSteps))
    params.set("mc steps", (sim.acceptSteps + sim.rejectSteps) / (sim.L*sim.L*sim.L))
  }

  def clear() {
    grid3d.clear()
  }

  def run() {
    val L = params.iget("L")
    sim = new HeisenbergNNSim(L)
    setSimParams()
    
    while (true) {
      Job.animate()
      
      for (i <- 0 until 1)
        sim.step()
    }
  }
}
