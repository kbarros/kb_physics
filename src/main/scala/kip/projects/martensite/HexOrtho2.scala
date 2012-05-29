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

class HexOrtho2Sim(params: Parameters) {
  
  val L = params.fget("L")
  val dx = params.fget("dx")
  val lp = (L / dx).toInt
  val lattice = new PeriodicLattice(lp, lp, 1, dx)
  import lattice._
  val random = new Random(params.iget("Random seed", 0))
  
  var dt = params.fget("dt")
  readParams()
  
  var energy = 0.0
  
  // x/y components of displacement vector u
  val ux = new Array[Double](N)
  val uy = new Array[Double](N)
  
  // time derivative of u
  val uxdt = new Array[Double](N)
  val uydt = new Array[Double](N)
  
  // symmetrized strains
  val e1 = new Array[Double](N)
  val e2 = new Array[Double](N)
  val e3 = new Array[Double](N)

  // stress components
  val sxx = new Array[Double](N)
  val syy = new Array[Double](N)
  val sxy = new Array[Double](N)
  
  val tau = -1d
  val a_b = 60d
  val a_s = 120d
  val K0 = 1d
  val eta0 = 0.1
  
  init()
  
  def readParams() {
    dt = params.fget("dt")
  }

  def init() {
    for (i <- 0 until N) {
      ux(i) = 0.5*random.nextGaussian()
      uy(i) = 0.5*random.nextGaussian()
      uxdt(i) = 0
      uydt(i) = 0
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
      val uyx = dX(uy, i)
      val uyy = dY(uy, i)

      // Symmetry adapted strains
      e1(i) = (uxx+uyy)/2
      e2(i) = (uxx-uyy)/2
      e3(i) = (uxy+uyx)/2
    }
    
    for (i <- 0 until N) {
      val lp_E2 = laplacian(e2, i)
      val lp_E3 = laplacian(e3, i)
      val E1 = e1(i) 
      val E2 = e2(i)
      val E3 = e3(i)
      
      //
      // F = + tau (e2^2 + e3^2) + 2 (e2^3 - 3 e2 e3^2) + (e2^2 + e3^2)^2
      //     + (a_b/2) e1^2
      //     + (K0/2) ( (grad e2)^2 + (grad e3)^2 )
      //
      
      energy += (dx*dx) * (
          tau * (sqr(E2) + sqr(E3)) +
          2 * (cube(E2) - 3*E2*sqr(E3)) +
          sqr(sqr(E2) + sqr(E3)) +
          (a_b/2) * sqr(E1) +
          (K0/2) * (grad2(e2, i) + grad2(e3, i))
      )
      
      // Gi = dF / dei
      val G1 = a_b*E1
      val G2 = 2*tau*E2 + 6*(sqr(E2)-sqr(E3)) + 4*(sqr(E2) + sqr(E3))*E2 - K0*lp_E2
      val G3 = 2*tau*E3 - 12*E2*E3            + 4*(sqr(E2) + sqr(E3))*E3 - K0*lp_E3
      
      sxx(i) = 0.5 * (G1 + G2)
      syy(i) = 0.5 * (G1 - G2)
      sxy(i) = 0.5 * (G3)
    }
    
    for (i <- 0 until N) {
      val sxxx = dX(sxx, i) 
      val syyy = dY(syy, i)
      val sxyx = dX(sxy, i)
      val sxyy = dY(sxy, i)
      
      // update velocity field
      uxdt(i) += dt*(sxxx+sxyy+eta0*laplacian(uxdt, i))
      uydt(i) += dt*(sxyx+syyy+eta0*laplacian(uydt, i))
      
      // update displacements
      ux(i) += dt*uxdt(i)
      uy(i) += dt*uydt(i)

//      println("del " + dt*uxdt(i))
    }
  }
}


object HexOrtho2 extends App {
  new Control(new HexOrtho2(), "Cubic to Tetragonal")  
}

class HexOrtho2 extends Simulation {
  val grid = new gfx.GridView()
  gfx.Utilities.frame(grid.canvas, w=300, h=300, title="e3")
  
  var sim: HexOrtho2Sim = _

  def load(c: Control) {
    params.add("L", 40.0)
    params.add("dx", 0.5)
    params.add("Random seed", 0)

    params.addm("dt", 0.05)
    
    params.add("energy")
  }

  def animate() {
    sim.readParams()
    val data = {
      new gfx.GridView.ArrayData(sim.lp, sim.lp, sim.e3, gfx.ColorGradient.blueRed(lo= -1.5, hi=1.5))
    }
    grid.display(data)
    params.set("energy", sim.energy)
//    println("max " + sim.e3.max)
  }

  def clear() {
  }

  def run() {
    sim = new HexOrtho2Sim(params);   

    while (true) {
      for (i <- 0 until 1)
        sim.step();
      Job.animate();
    }
  }
}
