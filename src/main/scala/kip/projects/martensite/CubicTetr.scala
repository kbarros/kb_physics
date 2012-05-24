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
  val uxdt = new Array[Double](N)
  val uydt = new Array[Double](N)
  val uzdt = new Array[Double](N)
  
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
  val eta0 = 0.2d
  
  init()
  
  def readParams() {
    dt = params.fget("dt")
  }

  def init() {
    for (i <- 0 until N) {
      ux(i) = 0.5*random.nextGaussian()
      uy(i) = 0.5*random.nextGaussian()
      uz(i) = 0.5*random.nextGaussian()
      uxdt(i) = 0
      uydt(i) = 0
      uzdt(i) = 0
    }
  }
  
  def step() {
    def sqr(x: Double) = x*x
    def cube(x: Double) = x*x*x
    def grad2(a: Array[Double], i: Int) = sqr(dX(a, i)) + sqr(dY(a, i)) + sqr(dZ(a, i))

    energy = 0
    
    for (i <- 0 until N) {
      val uxx = dX(ux, i)
      val uxy = dY(ux, i)
      val uxz = dZ(ux, i)
      val uyx = dX(uy, i)
      val uyy = dY(uy, i)
      val uyz = dZ(uy, i)
      val uzx = dX(uz, i)
      val uzy = dY(uz, i)
      val uzz = dZ(uz, i)

      // Symmetry adapted strains
      e1(i) = (uxx+uyy+uzz)/sqrt(3)
      e2(i) = (uxx-uyy)/(sqrt(2))
      e3(i) = (uxx+uyy-2*uzz)/(sqrt(6))
      e6(i) = (uxy+uyx)/(2*dx)
      e5(i) = (uxz+uzx)/(2*dx)
      e4(i) = (uyz+uzy)/(2*dx)
    }
    
    for (i <- 0 until N) {
      val lpe2 = laplacian(e2, i)
      val lpe3 = laplacian(e3, i)
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
      val G2 = 2*tau*E2 + 6*(sqr(E2)-sqr(E3)) + 4*(sqr(E2) + sqr(E3))*E2 - K0*lpe2
      val G3 = 2*tau*E3 - 12*E2*E3            + 4*(sqr(E2) + sqr(E3))*E3 - K0*lpe3
      
      sxx(i) = (G1/sqrt(3.)) + (G2/sqrt(2.)) + (G3/sqrt(6.))
      syy(i) = (G1/sqrt(3.)) - (G2/sqrt(2.)) + (G3/sqrt(6.))
      szz(i) = (G1/sqrt(3.)) - (2*G3/sqrt(6.))
      sxy(i) = a_s*E6
      sxz(i) = a_s*E5
      syz(i) = a_s*E4
    }
    
    for (i <- 0 until N) {
      val sxxx = dX(sxx, i) 
      val syyy = dY(syy, i)
      val szzz = dZ(szz, i)
      val sxyx = dX(sxy, i)
      val sxyy = dY(sxy, i)
      val sxzx = dX(sxz, i)
      val sxzz = dZ(sxz, i)
      val syzz = dZ(syz, i)
      val syzy = dY(syz, i)
      
      // divergence of sigma', for the drag
      val lpx = laplacian(uxdt, i)
      val lpy = laplacian(uydt, i)
      val lpz = laplacian(uzdt, i)
      
      // update velocity field
      uxdt(i) += dt*(sxxx+sxyy+sxzz+eta0*lpx)
      uydt(i) += dt*(sxyx+syyy+syzz+eta0*lpy)
      uzdt(i) += dt*(sxzx+syzy+szzz+eta0*lpz)
      
      // update displacements
      ux(i) += dt*uxdt(i)
      uy(i) += dt*uydt(i)
      uz(i) += dt*uzdt(i)
      
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
