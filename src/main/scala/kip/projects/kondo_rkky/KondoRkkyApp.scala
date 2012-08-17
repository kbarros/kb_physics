package kip.projects.kondo_rkky

import scikit.graphics.dim2.Grid
import scikit.jobs.Control
import scikit.jobs.Job
import scikit.jobs.Simulation
import scikit.jobs.params.DoubleValue
import scikit.jobs.params.Parameters
import scikit.numerics.Math2;
import scala.util.Random

import math._


class Rkky(L: Int, mu: Double, beta: Double) {
  val N = L*L
  
  val chi = new Array[Double](N)

  val eps = Array.tabulate[Double](N) { i =>
    val kx = i % L
    val ky = i / L
    -2 * (cos(2*Pi*kx/L) + cos(2*Pi*ky/L) + cos(2*Pi*(kx+ky)/L))
  }
  
  println("max %f max %f".format(eps.min, eps.max))
  
  val f = eps.map { e => 1 / (exp(beta * (e - mu)) + 1) }
  
  for (i_q <- 0 until N;
       i_k <- 1 until N) {
    
    val qx = i_q % L
    val qy = i_q / L
    
    val kx = i_k % L
    val ky = i_k / L

    val qkx = (qx+kx) % L
    val qky = (qy+ky) % L
    
    val i_qk = qky*L + qkx
    
    val de = eps(i_qk) - eps(i_k)
    
    def sqr(x: Double) = x*x
    
    val r1 = (beta/2) * (sqr(f(i_qk)) + sqr(f(i_k)) - f(i_qk) - f(i_k))
    val r2 = (f(i_qk) - f(i_k)) / de
    
//    if (abs(de) > 1e-6 && abs(de) < 1e-2)
//      println("de %g fs %g %g rs %g %g".format(de, f(i_qk), f(i_k), r1, r2))
      
    chi(i_q) += 3 * (if (abs(de) < 1e-5) r1 else r2)
  }
  chi(0) = 0
}


object KondoRkkyApp extends App {
  new Control(new KondoRkkyApp(), "Kondo Rkky Model");  
}

class KondoRkkyApp extends Simulation {
  val grid = new Grid("Susceptibility")

  def load(c: Control) {
//    grid.setScale(-1, +1)
    c.frame(grid)
    params.add("L", 50)
    params.addm("mu", -3.0)
    params.addm("beta", 100.0)
  }

  def animate() {
    val L = params.iget("L")
    val rkky = new Rkky(L, params.fget("mu"), params.fget("beta"))
    grid.registerData(L, L, rkky.chi)
  }

  def clear() {
    grid.clear()
  }

  def run() {
//    sim = new Annni2dSim(params);    

    while (true) {
//      for (i <- 0 until 1000)
//        sim.step();
      Job.animate();
    }
  }
}
