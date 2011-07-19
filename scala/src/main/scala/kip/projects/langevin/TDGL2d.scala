package kip.projects.langevin

import java.util.Random
import scikit.jobs.params.Parameters
//import scikit.numerics.fft.FFT2D
import kip.math.fft.FFTReal
import scikit.numerics.fn.Function2D

object TDGL2d {
  // constants from Cheng & Rutenberg
  val a1 = 3.0
  val a2 = 0.0
  
  // approximate energy at a single interface per unit area
  val surfaceEnergyDensity = 1.0 // 0.89
}


class TDGL2d(params: Parameters) {
  import TDGL2d._
  
  val random = new Random(params.iget("Random seed", 0))
  val L = params.fget("L")
  var t = 0.0
  var dt = params.fget("dt")
  val (lp, dx) = {
    var dx = params.fget("dx")
    val lp = (L / dx).toInt
    dx = L / lp
    params.set("dx", dx)
    (lp, dx)
  }
  val phi = new Array[Double](lp*lp)
  val scratch1 = new Array[Double](lp*lp)
  val scratch2 = new Array[Double](lp*lp)
  val fft = new FFTReal(Array(lp, lp), Some(Array(L, L)))
  
  randomize()
  
  
  def randomize() {
    for (i <- 0 until lp*lp)
      phi(i) = 0.1*random.nextGaussian()
    t = 0
  }
  
  def randomizeAndShift() {
    randomize()
    var sum = 0.0
    for (i <- 0 until lp*lp)
      sum += phi(i)
    for (i <- 0 until lp*lp)
      phi(i) -= sum / (lp*lp)
    
    sum = 0
    for (i <- 0 until lp*lp)
      sum += phi(i)
    if (sum > 1e-8 || sum < -1e-8) {
      System.out.println("Error in shifting!")
    }

    t = 0
  }

  def randomizeIsing() {
    val lp2 = lp*lp 
    val indices = new Array[Int](lp2)
    for (i <- 0 until lp2)
      indices(i) = i
    // randomize indices
    for (i <- 0 until lp2) {
      val r = i+random.nextInt(lp2-i)
      // swap indices(i) and indices(r)
      val that = indices(r)
      indices(r) = indices(i)
      indices(i) = that
    }
    for (i <- 0 until lp2) {
      phi(indices(i)) = if (i < lp2/2) -1 else +1
    }

    var sum = 0.0
    for (i <- 0 until lp2)
      sum += phi(i)
    if (sum != 0.0) {
      System.out.println("Error in randomizing!")
    }
    
    t = 0
  }
  
  def readParams(params: Parameters) {
    dt = params.fget("dt")
  }
  
  def r2k2(k: Array[Double]): Double = {
    val d = k.size
    require(d == 2)
    
    var k2 = 0.0
    for (i <- 0 until d)
      k2 += k(i)*k(i)

    var c = 1.0
    for (i <- 0 until d)
      c  *= d*k(i)*k(i)/k2

    val a = 0.5
    if (k2 == 0) 0 else ((1-a)+a*c)*k2
  }

  def simulate() {
    // scratch1 stores term proportional to phi
    fft.convolveWithFn(phi, scratch1) { k: Array[Double] =>
      val k2 = r2k2(k)
      val re = (1 + a1*dt + a2*dt*(-k2)) / (1 + (a1-1)*dt + (a2-1)*dt*(-k2))
      (re, 0)
    }
    
    // scratch2 stores term proportional to phi^3
    for (i <- 0 until lp*lp) {
      scratch2(i) = phi(i)*phi(i)*phi(i)
    }
    fft.convolveWithFn(scratch2, scratch2) { k: Array[Double] =>
      val k2 = r2k2(k)
      val re = dt / (1 + (a1-1)*dt + (a2-1)*dt*(-k2))
      (re, 0)
    }

    for (i <- 0 until lp*lp) {
      phi(i) = scratch1(i) - scratch2(i)
    }
    
    t += dt
  }
  
  // F = \int dx^2 [ - phi (1 + \del^2) phi / 2 + phi^4 / 4 ]  
  def freeEnergy: Double = {
    // scratch1 stores (1 + \del^2) phi
    fft.convolveWithFn(phi, scratch1) { k: Array[Double] =>
      val k2 = r2k2(k)
      (1-k2, 0)
    }
    var ret = 0.0
    for (i <- 0 until lp*lp) {
      ret += - phi(i)*scratch1(i)/2      // - phi (1 + \del^2) phi / 2 
      ret += phi(i)*phi(i)*phi(i)*phi(i) / 4	// phi^4 / 4
    }
    
    // shift energy so that it is zero in the ground state
    ret += lp*lp/4 
    
    ret*dx*dx
  }

}
