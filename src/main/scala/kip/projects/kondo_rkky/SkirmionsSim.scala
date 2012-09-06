package kip.projects.kondo_rkky

import scala.util.Random
import scala.math._
import kip.math.fft.FFTReal
import kip.math.Math._


class SkirmionsSim(val d: Int, val L: Int, var T: Double, var H: Double, var anisotropy: Double, var dt: Double) {
  val rand = new Random(System.currentTimeMillis())

  val dimensions = Array.fill(d)(L)
  val N = dimensions.product
  
  val sx = new Array[Double](N)
  val sy = new Array[Double](N)
  val sz = new Array[Double](N)
  randomizeSpins()
  
  val sbar = new Array[Double](N)

  val fft = new FFTReal(dim=dimensions)
  
  val kernel = {
    val ret = fft.allocFourierArray()
    for (i <- 0 until ret.size/2) {
      val k = fft.fourierVector(i)
      val k2 = k.map(sqr(_)).sum
      
      val w0 = L/8.0
      val q0 = 2*Pi*w0/L

      val inter =  sqr(k2 - sqr(q0)) 
      val inter0 = sqr(sqr(2*q0) - sqr(q0)) 

      ret(2*i+0) = min(inter, inter0)
      ret(2*i+1) = 0
    }
    ret
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
  
  // sum_q chi_q |S(q)|^2
  def energy(): Double = {
    var ret = 0.0
    for (si <- Seq(sx, sy, sz)) {
      fft.convolveWithRecip(si, kernel, sbar)
      for (i <- 0 until N) {
        ret += 0.5 * si(i) * sbar(i)
      }
    }
    ret / N
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
  }
}
