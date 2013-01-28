package kip.projects.kondo_rkky

import kip.graphics.Bounds3d
import kip.graphics.RetainedScene
import kip.math.Vec3
import scikit.graphics.dim2.Grid
import scikit.jobs.Control
import scikit.jobs.Job
import scikit.jobs.Simulation
import scikit.jobs.params.DoubleValue
import scikit.graphics.dim3.Grid3D


object Skirmions3dApp extends App {
  new Control(new Skirmions3dApp(), "3d Skyrmions")
}

class Skirmions3dApp extends Simulation {
  val grid3d = new Grid3D("Grid")
  val gridHH = new Grid3D("Hedgehog")

  var sim: SkirmionsSim = _
  val rs = new RetainedScene(Bounds3d(Vec3(0,0,0),Vec3(1,1,0)), sizew=600, sizeh=600, cameraDistance=1.2)

  def load(c: Control) {
    c.frame(grid3d, gridHH)
    params.add("L", 30)
    params.add("Len", 30.0)
    params.addm("T", 0.0)
    params.addm("H", 0.2)
    params.addm("anisotropy", 0.2)
    params.addm("dt", 0.1)
    params.addm("slice", new DoubleValue(0, 0, 1).withSlider)
    params.add("energy")
  }

  def animate() {
    sim.T = params.fget("T")
    sim.H = params.fget("H")
    sim.anisotropy = params.fget("anisotropy")
    sim.dt = params.fget("dt")
    
    val slice = math.min((params.fget("slice") * sim.L).toInt, sim.L-1)
    
    val field = new Array[Double](3*sim.L*sim.L)
    spinXYCrossSection(field, slice)
    val viz = new SpinViz(sim.L, sim.L)
    viz.drawSpins(field, rs)
    
    grid3d.setScale(-1, 1)
    val lp = sim.L
    grid3d.registerData(lp, lp, lp, sim.sz.map(- _))
    
    gridHH.setScale(-1, 1)
    val charge = Array.tabulate[Double](sim.L*sim.L*sim.L) { windingCharge(_) }
    gridHH.registerData(lp, lp, lp, charge)

    params.set("energy", sim.energy())
  }

  def clear() {
    grid3d.clear()
    gridHH.clear()
  }

  def spinXYCrossSection(field: Array[Double], z: Int) {
    for (i <- 0 until sim.L*sim.L) {
      val ip = i + z*sim.L*sim.L
      field(0 + 3*i) = sim.sx(ip)
      field(1 + 3*i) = sim.sy(ip)
      field(2 + 3*i) = sim.sz(ip)
    }
  }
  
  // spins in yz plane at fixed x index
  def spinYZCrossSection(field: Array[Double], x: Int) {
    var cnt = 0
    val L = sim.L
    for (z <- 0 until L;
         y <- 0 until L) {
      val i = z*L*L + y*L + x
      field(0 + 3*cnt) = sim.sy(i)
      field(1 + 3*cnt) = sim.sz(i)
      field(2 + 3*cnt) = sim.sx(i)
      cnt += 1
    }
  }
  
  def donutSpins() {
    val L = sim.L
    for (x <- 0 until L;
         y <- 0 until L;
         z <- 0 until L) {
      val i = z*L*L + y*L + x
      sim.sx(i) = 0
      sim.sy(i) = 0
      
      import kip.math.Math.sqr
      val r = math.sqrt(sqr(x-L/2) + sqr(y-L/2))
      if (z > L/8 && z < 2*L/8 && r > L/8 && r < 2*L/8)
        sim.sz(i) = -1
      else
        sim.sz(i) = 1
    }
  }
  
  def ballSpins(sigma: Double) {
    val L = sim.L
    for (x <- 0 until L;
         y <- 0 until L;
         z <- 0 until L) {
      val i = z*L*L + y*L + x
      
      val del = Vec3(x-L/2+0.5, y-L/2+0.5, z-L/2+0.5) * sim.dx
      val r = del.norm
      
      import kip.math.Math._
      val s_z = 1 - 2*math.exp(- r*r / (2 * sqr(sigma)))      
      
      val e_rho = Vec3(del.x, del.y, 0).normalize
      val s_rho = math.sqrt(1 - s_z*s_z)
      
      sim.sx(i) = e_rho.x * s_rho
      sim.sy(i) = e_rho.y * s_rho
      sim.sz(i) = s_z
    }
    
    sim.normalizeSpins()
  }
  
  def wrap(xp: Int) = (xp+sim.L)%sim.L
  
  def idx(xp: Int, yp: Int, zp: Int) = {
    wrap(zp)*sim.L*sim.L + wrap(yp)*sim.L + wrap(xp)
  }
  
  def spin(i: Int) = Vec3(sim.sx(i), sim.sy(i), sim.sz(i))
  
  def coords(i: Int) = {
    val L = sim.L
    val x = i % L
    val y = (i / L) % L
    val z = i / (L * L)
    (x, y, z)
  }
  
  // \int dx dy S dot (dS/dx) cross (dS/dy)
  def windingCharge(i: Int) = {
    val (x, y, z) = coords(i)
    
    def facexy(x: Int, y: Int, z: Int) = {
      val S1 = spin(idx(x, y, z))
      val S2 = spin(idx(x+1, y, z))
      val S3 = spin(idx(x, y+1, z))
      S1 dot (S2 cross S3)
    }
    
    def faceyz(x: Int, y: Int, z: Int) = {
      val S1 = spin(idx(x, y, z))
      val S2 = spin(idx(x, y+1, z))
      val S3 = spin(idx(x, y, z+1))
      S1 dot (S2 cross S3)
    }
    
    def facezx(x: Int, y: Int, z: Int) = {
      val S1 = spin(idx(x, y, z))
      val S2 = spin(idx(x, y, z+1))
      val S3 = spin(idx(x+1, y, z))
      S1 dot (S2 cross S3)
    }
    
    (facexy(x, y, z) - facexy(x, y, z+1) +
     faceyz(x, y, z) - faceyz(x+1, y, z) +
     facezx(x, y, z) - facezx(x, y+1, z)) / 8
  }

  def run() {
    val L = params.iget("L")
    val len = params.fget("Len")
    sim = new SkirmionsSim(d=3, L=L, len=len, T=params.fget("T"), H=params.fget("H"), anisotropy=params.fget("anisotropy"), dt=params.fget("dt"))
    // donutSpins()
//    ballSpins(2)

    val printEnergy = false
    if (printEnergy) {
      var sigma = 0.0
      for (i <- 0 until 120) {
        ballSpins(sigma)
        printf("%f %f\n", sigma, sim.energy())
        sigma += 0.05
      }
    }
    
//    val printPairEnergy = true
//    if (printPairEnergy) {
//      val sigma = 3.75
//      var disp = 1.0
//      for (i <- 0 until 120) {
//        ballSpins(sigma, Vec3)
//        printf("%f %f\n", sigma, sim.energy())
//        disp += 0.05
//      }
//    }
    
    while (true) {
      Job.animate()
      
//      for (i <- 0 until 20)
      sim.implicitStep()
//       sim.step()
    }
  }
}
