package kip.projects.kondo_rkky

import scala.util.Random
import scala.math._
import kip.math.fft.FFTReal
import kip.math.Math._


class SkirmionsSim(val d: Int, val L: Int, val len: Double, var T: Double, var H: Double, var anisotropy: Double, var dt: Double) {
  val dx = len / L
  val rand = new Random(System.currentTimeMillis())

  val dimensions = Array.fill(d)(L)
  val N = dimensions.product
  
  val sx = new Array[Double](N)
  val sy = new Array[Double](N)
  val sz = new Array[Double](N)
  randomizeSpins()
  
  val sbar = new Array[Double](N)

  val fft = new FFTReal(dim=dimensions, lenOption=Some(Array.fill(d)(len)))
  
  val kernel = {
    val ret = fft.allocFourierArray()
    for (i <- 0 until ret.size/2) {
      val k = fft.fourierVector(i)
      val k2 = k.map(sqr(_)).sum
      
      val ratio = 8 // lambda/L, ratio of skyrmion diameter to system length
      val q0 = 2*Pi / ratio
      val q1 = 2*q0
      val modulation = 1 / (1 + sqr(k2 / sqr(q1)))
      
//      val inter = sqr(sqrt(k2) - q0)

      ret(2*i+0) = sqr(k2 - sqr(q0)) * modulation
      ret(2*i+1) = 0
    }
    ret
  }
  
  val kernel2 = fft.allocFourierArray()
  
  def ferromagetizeSpins() {
    for (i <- 0 until N) {
      sx(i) = 0
      sy(i) = 0
      sz(i) = 1
    }
    normalizeSpins()
  }
  
  def randomizeSpins() {
    for (i <- 0 until N) {
      sx(i) = rand.nextGaussian()
      sy(i) = rand.nextGaussian()
      sz(i) = rand.nextGaussian()
    }
    normalizeSpins()
  }

  def normalizeSpins() {
    for (i <- 0 until N) {
      val len = math.sqrt(sx(i)*sx(i) + sy(i)*sy(i) + sz(i)*sz(i))
      sx(i) /= len
      sy(i) /= len
      sz(i) /= len
    }
  }
  
  // extensive energy
  def energy(): Double = {
    var ret = 0.0
    for (si <- Seq(sx, sy, sz)) {
      fft.convolveWithRecip(si, kernel, sbar)
      for (i <- 0 until N) {
        ret += 0.5 * si(i) * sbar(i)
      }
    }
    for (i <- 0 until N) {
      ret += - H * sz(i) - 0.5 * anisotropy * sz(i)*sz(i)
    }
    ret * dx*dx*dx
  }
  
  def implicitStep() {
    for (i <- 0 until kernel2.size/2) {
      kernel2(2*i+0) = 1 / (1 + dt * kernel(2*i+0))
      kernel2(2*i+1) = 0
    }

    for (si <- Seq(sx, sy, sz)) {
      fft.convolveWithRecip(si, kernel2, sbar)
      for (i <- 0 until N) {
        si(i) = sbar(i) + math.sqrt(dt * 2 * T) * rand.nextGaussian()
      }
    }
    for (i <- 0 until N) {
      sz(i) += dt * (H + anisotropy * sz(i))
    }
    normalizeSpins()
  }
  
  def step() {
    for (si <- Seq(sx, sy, sz)) {
      fft.convolveWithRecip(si, kernel, sbar)
      for (i <- 0 until N) {
        si(i) -= dt * sbar(i) + math.sqrt(dt * 2 * T) * rand.nextGaussian()
      }
    }
    for (i <- 0 until N) {
      sz(i) += dt * (H + anisotropy * sz(i))
    }
    
    normalizeSpins()
    
//    for (i <- 0 until N) {
//      if (i % L < L/2) {
//        sx(i) = 0
//        sy(i) = 0
//        sz(i) = 0
//      }
//    }
  }
}
