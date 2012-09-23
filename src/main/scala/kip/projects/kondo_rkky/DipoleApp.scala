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


object DipoleApp extends App {
  new Control(new DipoleApp(), "Dipole model")  
}

class DipoleApp extends Simulation {
  val grid3d = new Grid3D("S_x")
  var sim: DipoleSim2 = _
  val rs = new RetainedScene(Bounds3d(Vec3(0,0,0),Vec3(1,1,0)), sizew=600, sizeh=600, cameraDistance=1.2)
  
  def load(c: Control) {
    c.frame(grid3d)
    params.add("Lx", 30)
    params.add("Lz", 5)
    params.addm("T", 0.1)
    params.addm("H", 0.2)
    params.addm("J", 0.2)
    params.addm("anisotropy", 0.3)
    params.addm("dt", 0.5)
    params.addm("slice", new DoubleValue(0, 0, 1).withSlider)
    params.add("energy")
  }

  def animate() {
    sim.T = params.fget("T")
    sim.H = params.fget("H")
    sim.J = params.fget("J")
    sim.anisotropy = params.fget("anisotropy")
    sim.dt = params.fget("dt")
    
    val slice = math.min((params.fget("slice") * sim.Lz).toInt, sim.Lz-1)
    
    val field = new Array[Double](3*sim.Lx*sim.Lx)
    for (i <- 0 until sim.Lx*sim.Lx) {
      val ip = i + slice*sim.Lx*sim.Lx
      field(0 + 3*i) = sim.sx(ip)
      field(1 + 3*i) = sim.sy(ip)
      field(2 + 3*i) = sim.sz(ip)
    }
    val viz = new SpinViz(sim.Lx, sim.Lx)
    viz.drawSpins(field, rs)
    
    grid3d.setScale(-1, 1)
    grid3d.registerData(sim.Lx, sim.Lx, sim.Lz, sim.sz.map(- _))
    
    params.set("energy", sim.energy())
  }

  def clear() {
    grid3d.clear()
  }

  def run() {
    val Lx = params.iget("Lx")
    val Lz = params.iget("Lz")
//    sim = new DipoleSim(Lx=Lx, Ly=Lx, Lz=Lz, T=params.fget("T"), H=params.fget("H"), J=params.fget("J"), anisotropy=params.fget("anisotropy"), dt=params.fget("dt"))
    sim = new DipoleSim2(L=Lx, depth=Lz, T=params.fget("T"), H=params.fget("H"), J=params.fget("J"), anisotropy=params.fget("anisotropy"), dt=params.fget("dt"))
    
    while (true) {
      Job.animate()
      
      for (i <- 0 until 1)
        sim.step()
    }
  }
}
