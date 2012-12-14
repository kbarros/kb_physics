package kip.projects.kondo_rkky

import kip.graphics.Bounds3d
import kip.graphics.RetainedScene
import kip.math.Vec3
import scikit.graphics.dim2.Grid
import scikit.jobs.Control
import scikit.jobs.Job
import scikit.jobs.Simulation


object Skirmions2dApp extends App {
  new Control(new Skirmions2dApp(), "Skirmions")  
}

class Skirmions2dApp extends Simulation {
  val grid = new Grid("Spin")
  var sim: SkirmionsSim = _
  val rs = new RetainedScene(Bounds3d(Vec3(0,0,0),Vec3(1,1,0)), sizew=600, sizeh=600, cameraDistance=1.2)
  
  def load(c: Control) {
    c.frame(grid)
    params.add("L", 50)
    params.addm("T", 0.1)
    params.addm("H", 0.2)
    params.addm("anisotropy", 1.0)
    params.addm("dt", 0.01)
    params.add("energy")
    params.add("skyrmion charge")
  }
  
  def ballSpins(sigma: Double) {
    val L = sim.L
    for (x <- 0 until L;
         y <- 0 until L) {
      val i = y*L + x
      
      val del = Vec3(x-L/2+0.5, y-L/2+0.5, 0) * sim.dx
      val r = del.norm
      
      val s_z = 1 - 2*math.exp(- r*r / (2 * sigma))
      
      val s_rho = math.sqrt(1 - s_z*s_z)
      
      sim.sx(i) = del.x * s_rho / r
      sim.sy(i) = del.y * s_rho / r
      sim.sz(i) = s_z
    }
    
    sim.normalizeSpins()
  }


  def wrap(xp: Int) = (xp+sim.L)%sim.L
  
  def idx(xp: Int, yp: Int) = {
    wrap(yp)*sim.L + wrap(xp)
  }

  def coords(i: Int) = (i % sim.L, i / sim.L)
  
  def spin(i: Int) = Vec3(sim.sx(i), sim.sy(i), sim.sz(i))

    // \int dx dy S dot (dS/dx) cross (dS/dy)
  def windingCharge(i: Int) = {
    val (x, y) = coords(i)
    
    val S1 = spin(idx(x, y))
    val S2 = spin(idx(x+1, y))
    val S3 = spin(idx(x, y+1))
    S1 dot (S2 cross S3)
  }

  def animate() {
    sim.T = params.fget("T")
    sim.H = params.fget("H")
    sim.anisotropy = params.fget("anisotropy")
    sim.dt = params.fget("dt")
    
    val charge = Array.tabulate[Double](sim.L*sim.L) { windingCharge(_) }

//    grid.registerData(sim.L, sim.L, sim.sz.map(s => -s).toArray)
    grid.registerData(sim.L, sim.L, charge)
    
    val field = new Array[Double](3*sim.sx.size)
    for (i <- 0 until sim.sx.size) {
      field(0 + 3*i) = sim.sx(i)
      field(1 + 3*i) = sim.sy(i)
      field(2 + 3*i) = sim.sz(i)
    }
    val viz = new SpinViz(sim.L, sim.L)
    viz.drawSpins(field, rs)
    
    params.set("energy", sim.energy())
    params.set("skyrmion charge", charge.sum)
  }

  def clear() {
    grid.clear()
  }

  def run() {
    val L = params.iget("L")
    sim = new SkirmionsSim(d=2, L=L, len=L, T=params.fget("T"), H=params.fget("H"), anisotropy=params.fget("anisotropy"), dt=params.fget("dt"))
    ballSpins(3.75)
    
    while (true) {
      Job.animate()
      
      for (i <- 0 until 20)
        sim.implicitStep()
    }
  }
}
