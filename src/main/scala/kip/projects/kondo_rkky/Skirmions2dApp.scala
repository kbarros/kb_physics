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
    params.addm("T", 0.05)
    params.addm("H", 0.1)
    params.addm("anisotropy", 0.5)
    params.addm("dt", 0.2)
    params.add("energy")
  }

  def animate() {
    sim.T = params.fget("T")
    sim.H = params.fget("H")
    sim.anisotropy = params.fget("anisotropy")
    sim.dt = params.fget("dt")
    
    grid.registerData(sim.L, sim.L, sim.sz.map(s => -s).toArray)
    
    val field = new Array[Double](3*sim.sx.size)
    for (i <- 0 until sim.sx.size) {
      field(0 + 3*i) = sim.sx(i)
      field(1 + 3*i) = sim.sy(i)
      field(2 + 3*i) = sim.sz(i)
    }
    val viz = new SpinViz(sim.L, sim.L)
    viz.drawSpins(field, rs)
    
    params.set("energy", sim.energy())
  }

  def clear() {
    grid.clear()
  }

  def run() {
    val L = params.iget("L")
    sim = new SkirmionsSim(d=2, L=L, T=params.fget("T"), H=params.fget("H"), anisotropy=params.fget("anisotropy"), dt=params.fget("dt"))
    
    while (true) {
      Job.animate()
      
      for (i <- 0 until 20)
        sim.step()
    }
  }
}
