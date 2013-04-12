package kip.projects.martensite

import scikit.graphics.dim2.Grid
import scikit.jobs.Control
import scikit.jobs.Job
import scikit.jobs.Simulation
import scikit.jobs.params.DoubleValue
import scikit.jobs.params.Parameters
import scikit.numerics.Math2;
import scala.util.Random
import scala.math.sqrt
import kip.{graphics => gfx}

class CubicTetrSim(params: Parameters) {
  
  val L = params.fget("L")
  val dx = params.fget("dx")
  val lp = (L / dx).toInt
  val lattice = new PeriodicLattice(lp, lp, lp, dx)
  import lattice._
  val random = new Random(params.iget("Random seed", 0))
  
  var dt = params.fget("dt")
  readParams()
  
  var energy = 0.0
  
  // x/y/z components of displacement vector u
  val ux = new Array[Double](N)
  val uy = new Array[Double](N)
  val uz = new Array[Double](N)
  
  // time derivative of u
  val ux_t = new Array[Double](N)
  val uy_t = new Array[Double](N)
  val uz_t = new Array[Double](N)
  
  // symmetrized strains
  val e1 = new Array[Double](N)
  val e2 = new Array[Double](N)
  val e3 = new Array[Double](N)
  val e4 = new Array[Double](N)
  val e5 = new Array[Double](N)
  val e6 = new Array[Double](N)

  // stress components
  val sxx = new Array[Double](N)
  val syy = new Array[Double](N)
  val szz = new Array[Double](N)
  val sxy = new Array[Double](N)
  val sxz = new Array[Double](N)
  val syz = new Array[Double](N)
  
  val tau = -1d
  val a_b = 60d
  val a_s = 120d
  val K0 = 1d
  val eta0 = 10d // 0.2d
  
  init()
  
  def readParams() {
    dt = params.fget("dt")
  }

  def init() {
    for (i <- 0 until N) {
      ux(i) = 0.5*random.nextGaussian()
      uy(i) = 0.5*random.nextGaussian()
      uz(i) = 0.5*random.nextGaussian()
      ux_t(i) = 0
      uy_t(i) = 0
      uz_t(i) = 0
    }
  }
  
  def step() {
    def sqr(x: Double) = x*x
    def cube(x: Double) = x*x*x
    def grad2(a: Array[Double], i: Int) = sqr(dX(a, i)) + sqr(dY(a, i)) + sqr(dZ(a, i))

    energy = 0
    
    for (i <- 0 until N) {
      val ux_x = dX(ux, i)
      val ux_y = dY(ux, i)
      val ux_z = dZ(ux, i)
      val uy_x = dX(uy, i)
      val uy_y = dY(uy, i)
      val uy_z = dZ(uy, i)
      val uz_x = dX(uz, i)
      val uz_y = dY(uz, i)
      val uz_z = dZ(uz, i)
      
      // Symmetry adapted strains
      e1(i) = (ux_x+uy_y+uz_z)/sqrt(3)
      e2(i) = (ux_x-uy_y)/(sqrt(2))
      e3(i) = (ux_x+uy_y-2*uz_z)/(sqrt(6))
      e6(i) = (ux_y+uy_x)/2
      e5(i) = (ux_z+uz_x)/2
      e4(i) = (uy_z+uz_y)/2
    }
    
    for (i <- 0 until N) {
      val lp_E2 = laplacian(e2, i)
      val lp_E3 = laplacian(e3, i)
      val E1 = e1(i) 
      val E2 = e2(i)
      val E3 = e3(i)
      val E4 = e4(i) 
      val E5 = e5(i)
      val E6 = e6(i)
      
      //
      // F = + tau (e2^2 + e3^2) + 2 (e2^3 - 3 e2 e3^2) + (e2^2 + e3^2)^2
      //     + (a_b/2) e1^2 + (a_s/2) (e4^2 + e5^2 + e6^2)
      //     + (K0/2) ( (grad e2)^2 + (grad e3)^2 )
      //
      
      energy += (dx*dx*dx) * (
          tau * (sqr(E2) + sqr(E3)) +
          2 * (cube(E2) - 3*E2*sqr(E3)) +
          sqr(sqr(E2) + sqr(E3)) +
          (a_b/2) * sqr(E1) +
          (a_s/2) * (sqr(E4) + sqr(E5) + sqr(E6)) +
          (K0/2) * (grad2(e2, i) + grad2(e3, i))
      )
      
      // Gi = dF / dei
      val G1 = a_b*E1
      val G2 = 2*tau*E2 + 6*(sqr(E2)-sqr(E3)) + 4*(sqr(E2) + sqr(E3))*E2 - K0*lp_E2
      val G3 = 2*tau*E3 - 12*E2*E3            + 4*(sqr(E2) + sqr(E3))*E3 - K0*lp_E3
      
      sxx(i) = (G1/sqrt(3.0)) + (G2/sqrt(2.0)) + (G3/sqrt(6.0))
      syy(i) = (G1/sqrt(3.0)) - (G2/sqrt(2.0)) + (G3/sqrt(6.0))
      szz(i) = (G1/sqrt(3.0)) - (2*G3/sqrt(6.0))
      sxy(i) = a_s*E6
      sxz(i) = a_s*E5
      syz(i) = a_s*E4
    }
    
    for (i <- 0 until N) {
      val sxx_x = dX(sxx, i) 
      val syy_y = dY(syy, i)
      val szz_z = dZ(szz, i)
      val sxy_x = dX(sxy, i)
      val sxy_y = dY(sxy, i)
      val sxz_x = dX(sxz, i)
      val sxz_z = dZ(sxz, i)
      val syz_z = dZ(syz, i)
      val syz_y = dY(syz, i)
      
      // divergence of sigma', for the drag
      val lp_ux_t = laplacian(ux_t, i)
      val lp_uy_t = laplacian(uy_t, i)
      val lp_uz_t = laplacian(uz_t, i)
      
      // update velocity field
      ux_t(i) += dt*(sxx_x+sxy_y+sxz_z+eta0*lp_ux_t)
      uy_t(i) += dt*(sxy_x+syy_y+syz_z+eta0*lp_uy_t)
      uz_t(i) += dt*(sxz_x+syz_y+szz_z+eta0*lp_uz_t)
      
      // update displacements
      ux(i) += dt*ux_t(i)
      uy(i) += dt*uy_t(i)
      uz(i) += dt*uz_t(i)
      
//      println("del " + dt*uxdt(i))
    }
  }
}


object CubicTetr extends App {
  new Control(new CubicTetr(), "Cubic to Tetragonal")  
}

class CubicTetr extends Simulation {
  val grid = new gfx.GridView()
  gfx.Utilities.frame(grid.canvas, w=300, h=300, title="e3")
  
  var sim: CubicTetrSim = _

  def load(c: Control) {
    params.add("L", 20.0)
    params.add("dx", 1)
    params.add("Random seed", 0)

    params.addm("dt", 0.02)
    params.addm("slice", 0)
    
    params.add("energy")
  }

  def animate() {
    sim.readParams()
    val data = {
      val i = params.iget("slice")
      val n = sim.lp*sim.lp
      val a = sim.e3.slice(i*n, (i+1)*n)
      new gfx.GridView.ArrayData(sim.lp, sim.lp, a, gfx.ColorGradient.blueRed(lo= -1.5, hi=1.5))
    }
    grid.display(data)
    params.set("energy", sim.energy)
//    println("max " + sim.e3.max)
  }

  def clear() {
  }

  def run() {
    sim = new CubicTetrSim(params);   

    while (true) {
      for (i <- 0 until 1)
        sim.step();
      Job.animate();
    }
  }
}
