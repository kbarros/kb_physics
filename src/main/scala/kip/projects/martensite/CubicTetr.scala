package kip.projects.martensite

import scikit.graphics.dim2.Grid
import scikit.jobs.Control
import scikit.jobs.Job
import scikit.jobs.Simulation
import scikit.jobs.params.DoubleValue
import scikit.jobs.params.Parameters
import scikit.numerics.Math2
import scala.util.Random
import scala.math.sqrt
import kip.{graphics => gfx}
import scikit.graphics.dim3.Grid3D
import java.awt.Color

class CubicTetrSim(params: Parameters) {
  
  val random = new Random(params.iget("Random seed", 0))
  val L = params.fget("L")
  val dx = params.fget("dx")
  val dt = params.fget("dt")
  val tau = params.fget("tau")
  val a_b = params.fget("a_b")
  val a_s = params.fget("a_s")
  val K0 = params.fget("K0")
  val eta0 = params.fget("eta0")

  val lp = (L / dx).toInt
  val lattice = new PeriodicLattice(lp, lp, lp, dx)
  import lattice._
  
  // displacement vector u
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
  
  // elastic + GL energy
  val energy = new Array[Double](N) 
  
  def sqr(x: Double) = x*x
  def cube(x: Double) = x*x*x
  def hypot(x: Double, y: Double) = math.sqrt(x*x + y*y)
  def hypot(x: Double, y: Double, z: Double) = math.sqrt(x*x + y*y + z*z)

  init()

  def init() {
    for (i <- 0 until N) {
      val (x,y,z) = lattice.index2coord(i)
      val r = hypot(x-lp/2, y-lp/2, z-lp/2) * dx
      val sigma = 2.0
      val amp = 0*math.exp(-r*r / (2*sigma*sigma))
 
      val c = 0.5
      ux(i) = c*random.nextGaussian() + amp
      uy(i) = c*random.nextGaussian()
      uz(i) = c*random.nextGaussian()
      ux_t(i) = 0
      uy_t(i) = 0
      uz_t(i) = 0
   }
  }
  
  def step() {
    def grad2(a: Array[Double], i: Int) = sqr(dX(a, i)) + sqr(dY(a, i)) + sqr(dZ(a, i))
    
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
      energy(i) = (dx*dx*dx) * (
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
  
  class MartensiteGrid3D(title: String) extends Grid3D(title) {
    override def getColor(x: Int, y: Int, z: Int): Color = {
      val Array(w, h, d) = super.getDimensions()
      if (x < 0 || x >= w || y < 0 || y >= h || z < 0 || z >= d)
        new Color(0,0,0,0)
      else {
        val i = w*h*z+w*y+x
        val theta = math.atan2(sim.e2(i), sim.e3(i)) + math.Pi
        val r = math.hypot(sim.e2(i), sim.e3(i))
        
        val color = java.awt.Color.HSBtoRGB((theta/(2*math.Pi)).toFloat, 1f, math.min(1, r.toFloat))
        
        // new Color(0,0,sim.e2(i).toInt*255)
        new Color(color)
      }
    }
  }

  val gridMart = new MartensiteGrid3D("variants")
  
  // nval grid = new gfx.GridView()
  // gfx.Utilities.frame(grid.canvas, w=300, h=300, title="e3")
  val gridE1 = new Grid("e1")
  val gridE2 = new Grid("e2")
  val gridE3 = new Grid("e3")
  val gridEnergy = new Grid("energy")
    
  var sim: CubicTetrSim = _

  def load(c: Control) {
    c.frameTogether("Elasticities", gridE1, gridE2, gridE3, gridEnergy)
    c.frame(gridMart)
    
    params.add("Random seed", 0)
    params.add("L", 20.0)
    params.add("dx", 1.0)
    params.add("dt", 0.02)
    params.add("tau", -1.0)
    params.add("a_b", 60d)
    params.add("a_s", 120d)
    params.add("K0", 1.0d)
    params.add("eta0", 0.5d)

    params.addm("slice", new DoubleValue(0, 0, 1).withSlider)
     
    params.add("energy density")
  }

  def animate() {
    val sliceIdx = math.min((params.fget("slice") * sim.lp).toInt, sim.lp-1)
    def slice(a: Array[Double]) = {
      a.slice(sliceIdx*sim.lp*sim.lp, (sliceIdx+1)*sim.lp*sim.lp)
    }

//    val data = {
//      val i = params.iget("slice")
//      val n = sim.lp*sim.lp
//      val a = sim.e3.slice(i*n, (i+1)*n)
//      new gfx.GridView.ArrayData(sim.lp, sim.lp, a, gfx.ColorGradient.blueRed(lo= -1.5, hi=1.5))
//    }
//    grid.display(data)
    gridE1.registerData(sim.lp, sim.lp, slice(sim.e1))
    gridE2.registerData(sim.lp, sim.lp, slice(sim.e2))
    gridE3.registerData(sim.lp, sim.lp, slice(sim.e3))
    gridEnergy.registerData(sim.lp, sim.lp, slice(sim.energy))
    
    gridMart.registerData(sim.lp, sim.lp, sim.lp, sim.e2)
    
    params.set("energy density", sim.energy.sum/(sim.L*sim.L*sim.L))
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
