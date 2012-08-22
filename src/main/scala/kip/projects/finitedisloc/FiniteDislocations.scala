package kip.projects.finitedisloc


import scikit.util.Utilities.format
import scikit.graphics.dim2.Grid
import scikit.jobs.Control
import scikit.jobs.Job
import scikit.jobs.Simulation
import scikit.jobs.params.Parameters

import kip.math.fft.FFTReal
import java.util.Random

import smatrix.Constructors.realDbl._ 


object FiniteDislocSim {
  /*
  import kip.projects.finitedisloc._
  import FiniteDislocSim._
  val m = mat2(1, -6, 3, 5)
  polarDecomp(m)
  */
  type Mat2 = smatrix.Dense[smatrix.Scalar.RealDbl]
  
  def mat2(a11: Double, a21: Double, a12: Double, a22: Double): Mat2 = {
    val ret = dense(2, 2)
    ret.copy(Array(a11, a21, a12, a22))
    ret
  }
  
  def polarDecomp(m: Mat2): Mat2 = {
    val a = m.toArray
    val b = new Array[Double](4)
    val det = a(0)*a(3) - a(2)*a(1) // reminder: column major
    if (det <= 0) {
      println("warning: determinant=%g".format(det))
    }
    val sdet = if (det >= 0) 1 else -1
    b(0) = m(0, 0) + sdet * a(3) 
    b(1) = m(1, 0) - sdet * a(2)
    b(2) = m(0, 1) - sdet * a(1)
    b(3) = m(1, 1) + sdet * a(0)
    val n = math.sqrt(b(0)*b(0) + b(1)*b(1))
    b(0) /= n
    b(1) /= n
    b(2) /= n
    b(3) /= n
    val ret = dense(2, 2)
    ret.copy(b)
    ret
  }
}


class FiniteDislocSim(params: Parameters) {
  import FiniteDislocSim._
  
  val random = new Random(params.iget("Random seed", 0))
  val L = params.fget("L")
  var t = 0.0
  var dt = 0.0
  
  val (lp, dx) = {
    var dx = params.fget("dx")
    val lp = (L / dx).toInt
    dx = L / lp
    params.set("dx", dx)
    (lp, dx)
  }
  
  val lp2 = lp*lp
  
  val fft = new FFTReal(Array(lp, lp), Some(Array(L, L)))

  val alpha31  = new Array[Double](lp2)
  val alpha32  = new Array[Double](lp2)
  val d2       = new Array[Double](lp2)
  val d3       = new Array[Double](lp2)

  // Useful constants
  val Id = mat2(1, 0, 0, 1)
  val W1 = mat2(1, 0, 0, 1)  / 2.0
  val W2 = mat2(1, 0, 0, -1) / 2.0
  val W3 = mat2(0, 1, 1, 0)  / 2.0
  val W4 = mat2(0, -1, 1, 0) / 2.0
  
  val c1 = math.sqrt(2)
  val c2 = math.sqrt(2)
  val c3 = 1.0
  val A1 = 1.0
  val A2 = 1.0
  val A3 = 1.0
  
  val gamma1 = fft.allocFourierArray()
  val gamma2 = fft.allocFourierArray()
  val gamma3 = fft.allocFourierArray()
  val gamma4 = fft.allocFourierArray()
  for (i <- 0 until gamma1.size/2) {
    val kx: Double = fft.fourierVector(i)(0)
    val ky: Double = fft.fourierVector(i)(1)
    val k2 = kx*kx + ky*ky
    if (k2 > 0) {
      gamma1(2*i+0) = (kx*kx-ky*ky)/k2
      gamma2(2*i+0) = (-2*kx*ky)/k2
      gamma3(2*i+1) = ky/k2
      gamma4(2*i+1) = kx/k2
    }
  }
  
  // Temporary data
  val d1       = new Array[Double](lp2)
  val d4       = new Array[Double](lp2)
  val d1k      = fft.allocFourierArray()
  val d2k      = fft.allocFourierArray()
  val d3k      = fft.allocFourierArray()
  val d4k      = fft.allocFourierArray()
  val alpha31k = fft.allocFourierArray()
  val alpha32k = fft.allocFourierArray()
  val F = Array.fill[Mat2](lp2)(dense(2,2))
  val E = Array.fill[Mat2](lp2)(dense(2,2))
  val theta    = new Array[Double](lp2)

  val detG = new Array[Double](lp2)
  val e1   = new Array[Double](lp2)
  val e2   = new Array[Double](lp2)
  val e3   = new Array[Double](lp2)
  
  val energyDensity = new Array[Double](lp2)
  var energy: Double = _
  
  // partial derivatives
  val du_dd1 = new Array[Double](lp2)
  val du_dd2 = new Array[Double](lp2)
  val du_dd3 = new Array[Double](lp2)
  val du_dd4 = new Array[Double](lp2)
  
  // total derivatives
  val du_dd2_tot  = new Array[Double](lp2)
  val du_dd3_tot  = new Array[Double](lp2)
  val du_da31_tot = new Array[Double](lp2)
  val du_da32_tot = new Array[Double](lp2)
  
  
  def calculateFields() {
    fft.forwardTransform(d2, d2k)
    fft.forwardTransform(d3, d3k)
    fft.forwardTransform(alpha31, alpha31k)
    fft.forwardTransform(alpha32, alpha32k)
    
    // calculate d1k and d4k in fourier space
    for (i <- 0 until d1k.size/2) {
      val kx: Double = fft.fourierVector(i)(0)
      val ky: Double = fft.fourierVector(i)(1)
      val k2  = kx*kx + ky*ky
      if (k2 == 0) {
        d1k(2*i+0) = 0
        d1k(2*i+1) = 0
        d4k(2*i+0) = 0
        d4k(2*i+1) = 0
      }
      else {
        val k2m = kx*kx - ky*ky
        d1k(2*i+0) = (k2m*d2k(2*i+0) - 2*kx*ky*d3k(2*i+0) + ky*alpha31k(2*i+1) - kx*alpha32k(2*i+1)) / k2
        d1k(2*i+1) = (k2m*d2k(2*i+1) - 2*kx*ky*d3k(2*i+1) - ky*alpha31k(2*i+0) + kx*alpha32k(2*i+0)) / k2
        d4k(2*i+0) = (2*kx*ky*d2k(2*i+0) - k2m*d3k(2*i+0) - kx*alpha31k(2*i+1) - ky*alpha32k(2*i+1)) / k2
        d4k(2*i+1) = (2*kx*ky*d2k(2*i+1) - k2m*d3k(2*i+1) + kx*alpha31k(2*i+0) + ky*alpha32k(2*i+0)) / k2
      }
    }
    
    // calculate d1 and d4 in real space
    fft.backwardTransform(d1k, d1)
    fft.backwardTransform(d4k, d4)
    
    // build real space quantities F, E, {e_k}, energy, unconstrained energy derivatives
    energy = 0
    for (i <- 0 until lp2) {
      import kip.math.Math.sqr
      detG(i) = sqr(1-d1(i)) - sqr(d2(i)) - sqr(d3(i)) + sqr(d4(i))
      require(detG(i) >= 0, "Distortion is negative at %d, %g %g %g %g %g".format(i, detG(i), d1(i), d2(i), d3(i), d4(i)))
      require(detG(i) > 1e-6, "Distortion is singular "+detG(i))
      
      F(i)(0, 0) = (1 - d1(i) + d2(i)) / detG(i)
      F(i)(1, 0) = (d3(i) - d4(i)) / detG(i)
      F(i)(0, 1) = (d3(i) + d4(i)) / detG(i)
      F(i)(1, 1) = (1 - d1(i) - d2(i)) / detG(i)
      val FT = F(i).tran

      E(i) = (FT*F(i) - Id) / 2
      
      e1(i) = c1 * (E(i) * W1.tran).trace
      e2(i) = c2 * (E(i) * W2.tran).trace
      e3(i) = c3 * (E(i) * W3.tran).trace
      
      val Q = polarDecomp(F(i))
      theta(i) = -math.atan2(Q(1,0), Q(1,1))
      
      energyDensity(i) = (A1*sqr(e1(i)) + A2*sqr(e2(i)) + A3*sqr(e3(i))) * detG(i) / 2
      energy += energyDensity(i) * sqr(dx)
      
      // unconstrained derivatives
      val m1 = (FT*F(i)*W1*c1 - Id*e1(i)) * FT * (2*A1*e1(i)*detG(i))
      val m2 = (FT*F(i)*W2*c2 - Id*e2(i)) * FT * (2*A2*e2(i)*detG(i))
      val m3 = (FT*F(i)*W3*c3 - Id*e3(i)) * FT * (2*A3*e3(i)*detG(i))
      du_dd1(i) = ((m1+m2+m3)*W1.tran).trace
      du_dd2(i) = ((m1+m2+m3)*W2.tran).trace
      du_dd3(i) = ((m1+m2+m3)*W3.tran).trace
      du_dd4(i) = ((m1+m2+m3)*W4.tran).trace
    }
    
    // constrained total derivative
    fft.convolveWithRecip(du_dd1, gamma1, du_dd2_tot)
    fft.convolveWithRecip(du_dd1, gamma2, du_dd3_tot)
    fft.convolveWithRecip(du_dd1, gamma3, du_da31_tot)
    fft.convolveWithRecip(du_dd1, gamma4, du_da32_tot)
    for (i <- 0 until lp2) {
      du_dd2_tot(i) += du_dd2(i)
      du_dd3_tot(i) += du_dd3(i)
    }
  }
  
  def readParams(params: Parameters) {
    dt = params.fget("dt")
  }
  
  def randomize() {
    for (i <- 0 until lp2) {
      alpha31(i) = 0
      alpha32(i) = 0
      d2(i) = 0
      d3(i) = 0
    }
    
//    // simple dipole
//    val x = lp/2
//    val y = lp/2
//    alpha31(lp*y + x) = 1.0 / (dx*dx)
//    alpha31(lp*(y+3) + x) = -1.0 / (dx*dx)
    
    // dislocation wall
    val x1 = 3*lp/4-10
    val x2 = 1*lp/4+10
    for (y <- lp/5 until 4*lp/5 by 8) {
      alpha32(lp*(y+4) + x1) = +1.0 / (dx*dx)
      alpha32(lp*y + x2) = -1.0 / (dx*dx)
    }
  }
  
  def simulate() {
    calculateFields()
    
    for (i <- 0 until lp2) {
      d2(i) -= dt*du_dd2_tot(i)
      d3(i) -= dt*du_dd3_tot(i)
    }
    
    t += dt
  }
}

object FiniteDislocations extends App {
  new Control(new FiniteDislocations(), "Finite Deformation with Dislocations")
}

class FiniteDislocations extends Simulation {
  val grid1 = new Grid("Grid1")
  val grid2 = new Grid("Grid2")
  var sim: FiniteDislocSim = _

  def load(c: Control) {
    c.frame(grid1)
    c.frame(grid2)
    params.addm("dt", 0.1)
    params.add("dx", 1.0)
    params.add("L", 50.0)
    params.add("Random seed", 0)
    params.add("Time")
    params.add("Energy")
  }
  
  def animate() {
    sim.readParams(params)
    
//    grid1.setScale(-0.2, 0.2)
//    grid2.setScale(-0.2, 0.2)
    
    val lp = sim.lp
    grid1.registerData(lp, lp, sim.theta)
    grid2.registerData(lp, lp, sim.d4)
    
    println(sim.F(10))
    
    params.set("Time", format(sim.t))
    params.set("Energy", format(sim.energy))
  }
  
  def clear() {
    grid1.clear()
    grid2.clear()
  }
  
  def run() {
    sim = new FiniteDislocSim(params)
    sim.randomize()
    Job.animate()
    
    while (true) {
      for (i <- 0 until 1)
        sim.simulate()
      Job.animate()
    }
  }
}
