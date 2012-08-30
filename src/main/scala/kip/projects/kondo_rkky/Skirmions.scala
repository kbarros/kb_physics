package kip.projects.kondo_rkky

import scikit.graphics.dim2.Grid
import scikit.jobs.Control
import scikit.jobs.Job
import scikit.jobs.Simulation
import scikit.jobs.params.DoubleValue
import scikit.jobs.params.Parameters
import scikit.numerics.Math2
import scala.util.Random
import math._
import kip.math.fft.{FFTReal, FFTComplex}
import kip.math.Math._
import java.lang.System


class Skirmions(val L: Int, var T: Double, var H: Double, var dt: Double) {
  val rand = new Random(System.currentTimeMillis())

  val N = L*L
    
  val sx = new Array[Double](N)
  val sy = new Array[Double](N)
  val sz = new Array[Double](N)
  randomizeSpins()
  
  val sbar = new Array[Double](N)

  val fft = new FFTReal(dim=Array(L, L))
  
  val kernel = {
    val ret = fft.allocFourierArray()
    for (i <- 0 until ret.size/2) {
      val k = fft.fourierVector(i)
      val k2 = sqr(k(0)) + sqr(k(1))
      
      val w0 = L/8.0
      val q0 = 2*Pi*w0/L

      val inter =  sqr(k2 - sqr(q0)) 
      val inter0 = sqr(sqr(2*q0) - sqr(q0)) 

      ret(2*i+0) = min(inter, inter0)
      ret(2*i+1) = 0
    }
    ret
  }
  
  val chi = {
    val ret = new Array[Double](N)
    val chi_k = fft.uncompressFourierArray(kernel)
    for (i <- 0 until N)
      ret(i) = chi_k(2*i+0)
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
  
  def ferromagnetizeSpins() {
    for (i <- 0 until N) {
      sx(i) = 1
      sy(i) = 0
      sz(i) = 0
    }
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
      sy(i) -= dt * H
    }
    
    normalizeSpins()
  }
}


object SkirmionsApp extends App {
  new Control(new SkirmionsApp(), "Skirmions");  
}

class SkirmionsApp extends Simulation {
  val grid1 = new Grid("Chi")
  val grid2 = new Grid("Spin")
  var sim: Skirmions = _
  
  def load(c: Control) {
    c.frame(grid1, grid2)
    params.add("L", 100)
    params.addm("T", 0.1)
    params.addm("H", 0.1)
    params.addm("dt", 0.2)
    params.add("energy")
  }

  def animate() {
    sim.T = params.fget("T")
    sim.H = params.fget("H")
    sim.dt = params.fget("dt")
    
    grid1.registerData(sim.L, sim.L, sim.chi)
    grid2.registerData(sim.L, sim.L, sim.sy)
    
    params.set("energy", sim.energy())
  }

  def clear() {
    grid1.clear()
    grid2.clear()
  }

  def run() {
    sim = new Skirmions(L=params.iget("L"), T=params.fget("T"), H=params.fget("H"), dt=params.fget("dt"))
    
    while (true) {
      Job.animate();
      sim.step()
    }
  }
}
