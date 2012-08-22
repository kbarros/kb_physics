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


class Rkky(val L: Int, val mu: Double, val beta: Double, val seed: Int = 0) {
  val N = L*L

  val dt = 0.01
    
  val rand = new Random(seed)
  
  val sx = new Array[Double](N)
  val sy = new Array[Double](N)
  val sz = new Array[Double](N)
  randomizeSpins()
  
  val sbar = new Array[Double](N)
  
  val chi = {
    val ret = new Array[Double](N)
    
    // tight binding eigenvalues
    val eps = Array.tabulate[Double](N) { i =>
      val kx = i % L
      val ky = i / L
      -2 * (cos(2*Pi*kx/L) + cos(2*Pi*ky/L) + cos(2*Pi*(kx+ky)/L))
    }

    // fermi function of tight binding eigenvalues
    val f = eps.map { e => 1 / (exp(beta * (e - mu)) + 1) }
    
    // susceptibility
    for (i_q <- 0 until N) {
      ret(i_q) = 0
      for (i_k <- 0 until N) {
        val qx = i_q % L
        val qy = i_q / L
        val kx = i_k % L
        val ky = i_k / L

        val qkx = (qx+kx) % L
        val qky = (qy+ky) % L

        val i_qk = qky*L + qkx

        val de = eps(i_qk) - eps(i_k)

        val r1 = (beta/2) * (sqr(f(i_qk)) + sqr(f(i_k)) - f(i_qk) - f(i_k))
        val r2 = (f(i_qk) - f(i_k)) / de

        ret(i_q) += (if (abs(de) < 1e-5) r1 else r2)
      }
    }
    
    ret
  }

  val fft = new FFTReal(dim=Array(L, L))
  
  val kernel = fft.allocFourierArray()
  for (i <- 0 until kernel.size/2) {
    val k = fft.fourierIndices(i)
    val kx = (k(0) + L) % L
    val ky = (k(1) + L) % L
    val i_k = ky*L + kx
    kernel(2*i+0) = chi(i_k)
    kernel(2*i+1) = 0
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
        ret += - 0.5 * si(i) * sbar(i)
      }
    }
    ret
  }
  
  def step() {
    for (si <- Seq(sx, sy, sz)) {
      fft.convolveWithRecip(si, kernel, sbar)
      for (i <- 0 until N) {
        si(i) -= dt * sbar(i) + math.sqrt(dt / beta) * rand.nextGaussian()
      }
    }
    normalizeSpins()
  }
}


object KondoRkkyApp extends App {
  new Control(new KondoRkkyApp(), "Kondo Rkky Model");  
}

class KondoRkkyApp extends Simulation {
  val grid1 = new Grid("Chi")
  val grid2 = new Grid("Spin")
  var sim: Rkky = _
  
  def load(c: Control) {
    c.frame(grid1, grid2)
    params.add("L", 10)
    params.add("mu", -3.0)
    params.add("beta", 100.0)
    params.add("energy")
  }

  def animate() {
    grid1.registerData(sim.L, sim.L, sim.chi)
//    grid.registerData(sim.L, sim.L, sim.sx)
    
    params.set("energy", sim.energy())
  }

  def clear() {
    grid1.clear()
    grid2.clear()
  }

  def run() {
    sim = new Rkky(params.iget("L"), params.fget("mu"), params.fget("beta"))
    
    println(sim.energy())
    
    while (true) {
      Job.animate();
      for (i <- 0 until 10)
        sim.step()
    }
  }
}
