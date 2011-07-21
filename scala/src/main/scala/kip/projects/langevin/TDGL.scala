package kip.projects.langevin

import java.util.Random

import kip.math.fft.FFTReal
import scikit.jobs.params.Parameters

object TDGL {
  // constants from Cheng & Rutenberg
  val a1 = 3.0
  val a2 = 0.0
  
  // approximate energy at a single interface per unit area
  val surfaceEnergyDensity = 1.0 // 0.89
}


class TDGL(params: Parameters, dimensions: Int) {
  import TDGL._
  
  val random = new Random(params.iget("Random seed", 0))
  val L = params.fget("L")
  var t = 0.0
  var dt = 0.0
  var r = 0.0
  
  val (lp, dx) = {
    var dx = params.fget("dx")
    val lp = (L / dx).toInt
    dx = L / lp
    params.set("dx", dx)
    (lp, dx)
  }
  
  val dim = Array.fill(dimensions)(lp)
  val len = Array.fill(dimensions)(L)
  val lpd = dim.product
  
  val phi = new Array[Double](lpd)
  val scratch1 = new Array[Double](lpd)
  val scratch2 = new Array[Double](lpd)
  val fft = new FFTReal(dim, Some(len))
  
  val kernel1 = fft.allocFourierArray()
  val kernel2 = fft.allocFourierArray()
  
  readParams(params)
  randomize()
  
  
  def randomize() {
    for (i <- 0 until lpd)
      phi(i) = 0.1*random.nextGaussian()
    t = 0
  }
  
  def randomizeAndShift() {
    randomize()
    var sum = phi.sum
    for (i <- 0 until lpd)
      phi(i) -= sum / (lpd)
    
    sum = phi.sum
    if (sum > 1e-8 || sum < -1e-8)
      println("Error in shifting!")
  }

  def randomizeIsing() {
    val indices = new Array[Int](lpd)
    for (i <- 0 until lpd)
      indices(i) = i
    // randomize indices
    for (i <- 0 until lpd) {
      val r = i+random.nextInt(lpd-i)
      // swap indices(i) and indices(r)
      val that = indices(r)
      indices(r) = indices(i)
      indices(i) = that
    }
    for (i <- 0 until lpd) {
      phi(indices(i)) = if (i < lpd/2) -1 else +1
    }

    var sum = 0.0
    for (i <- 0 until lpd)
      sum += phi(i)
    if (sum != 0.0) {
      println("Error in randomizing!")
    }
    
    t = 0
  }
  
  def readParams(params: Parameters) {
    dt = params.fget("dt")
    val rp  = params.fget("r")
    
    if (r != rp) {
      r = rp
      fft.tabulateFourierArray(kernel1) { k: Array[Double] =>
        val k2 = r2k2(k)
        val re = (1 + a1*dt + a2*dt*(-k2)) / (1 + (a1-1)*dt + (a2-1)*dt*(-k2))
        (re, 0)
      }
      fft.tabulateFourierArray(kernel2) { k: Array[Double] =>
        val k2 = r2k2(k)
        val re = dt / (1 + (a1-1)*dt + (a2-1)*dt*(-k2))
        (re, 0)
      }
    }
  }
  
  def r2k2(k: Array[Double]): Double = {
    require(k.size == dimensions)
    
    var k2 = 0.0
    for (i <- 0 until dimensions)
      k2 += k(i)*k(i)

    var c = 1.0
    for (i <- 0 until dimensions)
      c *= dimensions*k(i)*k(i)/k2

    val a = 0.5
    if (k2 == 0) 0 else r*r*((1-a)+a*c)*k2
  }
  
  def simulate() {
    // scratch1 stores term proportional to phi
    fft.convolveWithRecip(phi, scratch1)(kernel1)
    
    // scratch2 stores term proportional to phi^3
    for (i <- 0 until lpd) {
      scratch2(i) = phi(i)*phi(i)*phi(i)
    }
    fft.convolveWithRecip(scratch2, scratch2)(kernel2)

    for (i <- 0 until lpd) {
      phi(i) = scratch1(i) - scratch2(i)
    }
    
    t += dt
  }
  
  // F = \int dx^2 [ - phi (1 + \del^2) phi / 2 + phi^4 / 4 ]  
  def freeEnergy: Double = {
    // scratch1 stores (1 + \del^2) phi
    fft.convolveWithRecipFn(phi, scratch1) { k: Array[Double] =>
      val k2 = r2k2(k)
      (1-k2, 0)
    }
    var ret = 0.0
    for (i <- 0 until lpd) {
      ret += - phi(i)*scratch1(i)/2      // - phi (1 + \del^2) phi / 2 
      ret += phi(i)*phi(i)*phi(i)*phi(i) / 4	// phi^4 / 4
    }
    
    // shift energy so that it is zero in the ground state
    ret += lpd/4 
    
    ret*math.pow(dx, dimensions)
  }
}
