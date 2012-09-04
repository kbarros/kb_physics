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
import kip.math.Vec3
import kip.graphics.RetainedScene
import kip.graphics.Bounds3d
import javax.media.opengl.awt.GLJPanel
import javax.swing.SwingUtilities


class Skirmions(val L: Int, var T: Double, var H: Double, var anisotropy: Double, var dt: Double) {
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
      sx(i) = 0
      sy(i) = 0
      sz(i) = 1
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
      sz(i) += dt * (H + anisotropy * sz(i))
    }
    
    normalizeSpins()
  }
}


class SpinViz(w: Int, h: Int) {
  val latDel1 = Vec3(0.5, -0.5*math.sqrt(3), 0)
  val latDel2 = Vec3(0.5, +0.5*math.sqrt(3), 0)
  val latDel3 = latDel1 + latDel2

  val lat0 = Vec3(0, 0.5*math.sqrt(3)*(h-1), 0)
  val latN = lat0 + latDel1*(h-1) + latDel3*(w-1)
  
  val bds = Bounds3d(lat0, latN)

  def readSpin(x: Int, y: Int, field: Array[Double]): Vec3 = {
    require(x < w && y < w)
    val sx = field(0 + 3*(x + w*y))
    val sy = field(1 + 3*(x + w*y))
    val sz = field(2 + 3*(x + w*y))
    Vec3(sx, sy, sz)
  }
  
  def readPos(x: Int, y: Int): Vec3 = {
    lat0 + latDel3*x + latDel1*(h-1-y)
  }
  
  val spinPos: Array[Vec3] = {
    val ret = for (y <- (h-1) to 0 by -1;
         x <- 0 until w) yield {
      readPos(x, y)
    }
    ret.toArray
  }

  def spinDir(field: Array[Double]): Array[Vec3] = {
    val ret = for (y <- (h-1) to 0 by -1;
                   x <- 0 until w) yield {
      readSpin(x, y, field)
    }
    ret.toArray
  }

  def drawSpins(field: Array[Double], rs: RetainedScene) {

    val sd = spinDir(field)
    
    val arrows = for (i <- 0 until w*w) yield {
      val pos = spinPos(i) + Vec3(0, 0, 1)
      val spin = sd(i)
      val delta = spin*1.5
      val width = 0.3
      
      import java.awt.Color._
      val gray = new java.awt.Color(0, 0, 0, 50)

      new RetainedScene.Arrow(pos, delta, width, color1=ORANGE, color2=RED)
    }
    
    rs.bds = bds
    rs.drawables = Vector()
    rs.drawables ++= arrows
    
    SwingUtilities.invokeLater(new Runnable() {
      def run() = rs.display()
    })
  }

}

object SkirmionsApp extends App {
  val canvas = new GLJPanel() // new GLCanvas()
  
  new Control(new SkirmionsApp(), "Skirmions");  
}

class SkirmionsApp extends Simulation {
  val grid1 = new Grid("Chi")
  val grid2 = new Grid("Spin")
  var sim: Skirmions = _
  val rs = new RetainedScene(Bounds3d(Vec3(0,0,0),Vec3(1,1,0)), sizew=800, sizeh=600, cameraDistance=0.9)
  
  def load(c: Control) {
    c.frame(grid1, grid2)
    params.add("L", 50)
    params.addm("T", 0.05)
    params.addm("H", 0.1)
    params.addm("anisotropy", 0.5)
    params.addm("dt", 0.2)
    params.add("energy")
  }

  def animate() {
    sim.T = params.fget("T")
    sim.H = params.fget("H")
    sim.anisotropy = params.fget("anisotropy")
    sim.dt = params.fget("dt")
    
    grid1.registerData(sim.L, sim.L, sim.chi)
    grid2.registerData(sim.L, sim.L, sim.sz.map(s => -s).toArray)
    
    val field = new Array[Double](3*sim.sx.size)
    for (i <- 0 until sim.sx.size) {
      field(0 + 3*i) = sim.sx(i)
      field(1 + 3*i) = sim.sy(i)
      field(2 + 3*i) = sim.sz(i)
    }
    
    val viz = new SpinViz(sim.L, sim.L)
    viz.drawSpins(field, rs)
    
    params.set("energy", sim.energy())
  }

  def clear() {
    grid1.clear()
    grid2.clear()
  }

  def run() {
    val L = params.iget("L")
    sim = new Skirmions(L=L, T=params.fget("T"), H=params.fget("H"), anisotropy=params.fget("anisotropy"), dt=params.fget("dt"))
    
    while (true) {
      Job.animate()
      
      for (i <- 0 until 20)
        sim.step()
    }
  }
}
