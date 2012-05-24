package kip.projects.martensite

import scikit.graphics.dim2.Grid
import scikit.jobs.Control
import scikit.jobs.Job
import scikit.jobs.Simulation
import scikit.jobs.params.DoubleValue
import scikit.jobs.params.Parameters
import scikit.numerics.Math2;
import scala.util.Random
import kip.{graphics => gfx}


class HexOrthoSim(params: Parameters) {
  
  val L = params.iget("L")
  val dx = params.fget("dx")
  val lp = (L / dx).toInt 
  val lattice = new PeriodicLattice(lx=lp, ly=lp, lz=1, dx)
  import lattice._
  
  val random = new Random(params.iget("Random seed", 0))
  
  var zeta: Double = _
  var tau: Double = _
  var rho: Double = _
  var dt: Double = _
  
  var energy = 0.0
  
  // x/y components of displacement vector u
  val u1 = new Array[Double](N)
  val u2 = new Array[Double](N)
  
  // time derivative of u
  val udot_1 = new Array[Double](N)
  val udot_2 = new Array[Double](N)
  
  // symmetrized strains
  // e1 = (u_1,1 + u_2,2) / 2
  // e2 = (u_1,1 - u_2,2) / 2
  // e3 = (u_1,2 + u_2,1) / 2  
  val e1 = new Array[Double](N)
  val e2 = new Array[Double](N)
  val e3 = new Array[Double](N)
  
  // G1 = dF/de1 = zeta e1
  // G2 = dF/de2 = tau e2 - b e3^2  + c (e2^2 + e3^2) e2 - grad^2 e2
  // G3 = dF/de2 = tau e3 - b e2 e3 + c (e2^2 + e3^2) e3 - grad^2 e3
  val G1 = new Array[Double](N)
  val G2 = new Array[Double](N)
  val G3 = new Array[Double](N)

  readParams()
  init()

  def readParams() {
    zeta = params.fget("zeta")
    tau = params.fget("tau")
    rho = params.fget("rho")
    dt = params.fget("dt")
  }

  def init() {
    for (i <- 0 until N) {
      u1(i) = 0.1*random.nextGaussian()
      u2(i) = 0.1*random.nextGaussian()
      udot_1(i) = 0
      udot_2(i) = 0
    }
  }
  
  def step() {
    energy = 0
    
    for (i <- 0 until N) {
      // symmetrized strain fields
      e1(i) = (dX(u1, i) + dY(u2, i)) / 2
      e2(i) = (dX(u1, i) - dY(u2, i)) / 2
      e3(i) = (dY(u1, i) + dX(u2, i)) / 2
    }
    
    def sqr(x: Double) = x*x
    def cube(x: Double) = x*x*x
    def grad2(a: Array[Double], i: Int) = sqr(dX(a, i)) + sqr(dY(a, i))

    val r2 = 1.0
    
    for (i <- 0 until N) {
      // E = (zeta/2) e1^2 + (tau/2)(e2^2 + e3^2) - (b/3) (e2^3 - 3 e2 e3^2) + (c/4) (e2^2 + e3^2)^2 + (1/2) ((grad e_2)^2 + (grad e_3)^2)       
      val b = 3e3
      val c = 2e6
      energy += (dx*dx) * (
        (zeta/2d)*sqr(e1(i)) +
          (tau/2d)*(sqr(e2(i)) + sqr(e3(i))) +
        - (b/3d)*(cube(e2(i)) - 3*e2(i)*sqr(e3(i))) +
          (c/4d)*sqr(sqr(e2(i)) + sqr(e3(i))) +
          r2*(1/2d)*(grad2(e2, i) + grad2(e3, i))
      )
          
      // derivatives of free energy wrt e_i
      G1(i) = zeta * e1(i)
      G2(i) = tau*e2(i) - b*(sqr(e2(i))-sqr(e3(i))) + c*(sqr(e2(i))+sqr(e3(i)))*e2(i) - r2*laplacian(e2, i)
      G3(i) = tau*e3(i) + 2*b*e2(i)*e3(i)           + c*(sqr(e2(i))+sqr(e3(i)))*e3(i) - r2*laplacian(e3, i)
    }
    
//    for (i <- 0 until N) {
//      val ee = e2
//      energy += (dx*dx) * ( (0.25)*sqr(sqr(ee(i))-1) + (0.1*0.5)*(grad2(e1, i) + grad2(e2, i) + grad2(e3, i)) )
//      
//      // derivatives of free energy wrt e_i
//      G1(i) = - 0.1*laplacian(e1, i)
//      G2(i) = (sqr(ee(i))-1)*ee(i) - 0.1*laplacian(e2, i)
//      G3(i) = - 0.1*laplacian(e3, i)
//    }
    
    
    for (i <- 0 until N) {
      // potential driving force sigma_ij,j where sigma_ij = (dF/du_ij)
      val sigma_1jj = dX(G1, i) + dX(G2, i) + dY(G3, i)
      val sigma_2jj = dY(G1, i) - dY(G2, i) + dX(G3, i)
      
      // friction driving force sigma'_ij,j
      val A2p = 10.0
      val sigmap_1jj = A2p * laplacian(udot_1, i)/2
      val sigmap_2jj = A2p * laplacian(udot_2, i)/2
      
      // integrate rho (u_dot_dot) = div (sigma + sigma') 
      udot_1(i) += (dt / rho) * (sigma_1jj + sigmap_1jj)
      udot_2(i) += (dt / rho) * (sigma_2jj + sigmap_2jj)
      
      u1(i) += dt * udot_1(i)
      u2(i) += dt * udot_2(i)
    }
  }
}


object HexOrtho extends App {
  new Control(new HexOrtho(), "Hex Ortho")  
}

class HexOrtho extends Simulation {
  val grid1 = new gfx.GridView()
  val grid2 = new gfx.GridView()
  val grid3 = new gfx.GridView()
  gfx.Utilities.frame(grid1.canvas, w=300, h=300, title="e1")
  gfx.Utilities.frame(grid2.canvas, w=300, h=300, title="e2")
  gfx.Utilities.frame(grid3.canvas, w=300, h=300, title="e3")
  
  var sim: HexOrthoSim = _
  
  def load(c: Control) {
    params.add("L", 1)
    params.add("dx", 0.1)
    params.add("Random seed", 0)
    
    params.addm("zeta", new DoubleValue(0, 0, 2e3).withSlider())
    params.addm("tau",  new DoubleValue(-50, -50, 2).withSlider())
    params.addm("rho", 1.0)
    params.addm("dt", 1e-6)
    
    params.add("energy")
  }

  def animate() {
    sim.readParams()
    grid1.display(new gfx.GridView.ArrayData(sim.lp, sim.lp, sim.u1, gfx.ColorGradient.blueRed(lo= -0.05, hi=0.05)))
    grid2.display(new gfx.GridView.ArrayData(sim.lp, sim.lp, sim.u2, gfx.ColorGradient.blueRed(lo= -0.05, hi=0.05)))
    grid3.display(new gfx.GridView.ArrayData(sim.lp, sim.lp, sim.e3, gfx.ColorGradient.blueRed(lo= -0.05, hi=0.05)))
    params.set("energy", sim.energy)
  }

  def clear() {
    grid1.clear()
  }

  def run() {
    sim = new HexOrthoSim(params);    

    while (true) {
      for (i <- 0 until 1)
        sim.step();
      Job.animate();
    }
  }
}
