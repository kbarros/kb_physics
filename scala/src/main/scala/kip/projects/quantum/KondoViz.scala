package kip.projects.quantum

import net.liftweb.json
import java.io.File
import scikit.graphics.dim2.Grid
import kip.math.Vec3
import kip.enrich._
import kip.math.fft.FFTComplex
import kip.math.fft.FFTReal

object KondoViz extends App {
  val dir = args(0)
  val conf = KondoApp.readConfig(new File(dir+"/cfg.json"))
  val w = conf.w
  val h = conf.h
  val dumpdir = new java.io.File(dir+"/dump")
  require(dumpdir.isDirectory(), "Cannot load directory %s".format(dumpdir))
  
  def splitFieldComponents(field: Array[R]): Seq[Array[R]] = {
    val x = new Array[R](w*h)
    val y = new Array[R](w*h)
    val z = new Array[R](w*h)
    for (i <- 0 until w*h) {
      x(i) = field(3*i+0)
      y(i) = field(3*i+1)
      z(i) = field(3*i+2)
    }
    Seq(x, y, z)
  }
  
  def fftForward(a: Array[R]): Array[Double] = {
    val fft = new FFTReal(Array(h, w))
    val ap = fft.allocFourierArray()
    fft.forwardTransform(a.map(_.toDouble), ap)
    fft.uncompressFourierArray(ap)
  }
  
  def readSpin(x: Int, y: Int, field: Array[R]): Vec3 = {
    val sx = field(0 + x*3 + y*3*w)
    val sy = field(1 + x*3 + y*3*w)
    val sz = field(2 + x*3 + y*3*w)
    Vec3(sx, sy, sz)
  }

  import kip.graphics._
  import kip.math.Vec3
  val bds = Bounds3d(Vec3(0, 0, 0), Vec3(w-1, h-1, 0))
  val viz = new RetainedScene(bds)
  def drawSpins(field: Array[R]) {
    val arrows = for (y <- 0 until h;
                      x <- 0 until w) yield {
//                      if x % 2 == 0 && y % 2 == 1) yield {
      //val s = readSpin(x, y, field) * 0.8
      val s = readSpin(x, y, field) * 1.6
      val origin = Vec3(x, y, 0)
      val delta  = Vec3(s.x, s.y, s.z)
      new RetainedScene.Arrow(origin, delta, width=0.1, color1=java.awt.Color.RED, color2=java.awt.Color.BLUE)
    }
    viz.drawables = Vector(new RetainedScene.Cuboid(bds))
    viz.drawables ++= arrows
    viz.display()
  }

  val grid = new Grid("Order parameter")
  grid.setScale(-0.77, 0.77)
  scikit.util.Utilities.frame(grid.getComponent(), grid.getTitle())
  val gridData = new Array[Double]((2*w)*(2*h))
  def drawGrid(field: Array[R]) {
    for (y <- 0 until h;
         x <- 0 until w) yield {
      //
      // s4 - s3
      //  | / |
      // s1 - s2
      //
      val s1 = readSpin((x+0)%w, (y+0)%h, field)
      val s2 = readSpin((x+1)%w, (y+0)%h, field)
      val s3 = readSpin((x+1)%w, (y+1)%h, field)
      val s4 = readSpin((x+0)%w, (y+1)%h, field)
      gridData((2*y+0)*(2*w) + 2*x+0) = s1 dot (s2 cross s3)
      gridData((2*y+1)*(2*w) + 2*x+0) = s1 dot (s3 cross s4)
      gridData((2*y+0)*(2*w) + 2*x+1) = s1 dot (s2 cross s3)
      gridData((2*y+1)*(2*w) + 2*x+1) = s1 dot (s3 cross s4)
    }
    grid.registerData(2*w, 2*h, gridData)
  }
  
  val gridFft = new Grid("FFT")
  scikit.util.Utilities.frame(gridFft.getComponent(), gridFft.getTitle())
  val gridFftData = new Array[Double](w*h)
  def drawGridFft(field: Array[R]) {
    val Seq(sx, sy, sz) = splitFieldComponents(field)
    val fx = fftForward(sx)
    val fy = fftForward(sy)
    val fz = fftForward(sz)
    gridFftData.transform(_ => 0)
    for (i <- 0 until h*w) {
      import kip.math.Math.sqr
      for (f <- Seq(fx, fy, fz)) {
        gridFftData(i) += sqr(f(2*i+0)) + sqr(f(2*i+1))
      }
    }
    gridFft.registerData(w, h, gridFftData)
  }
  
//  val plot = KPM.mkPlot("Integrated rho")
//  def drawDensity(moments: Array[R]) {
//    val order = moments.size
//    val range = KPM.range(5*order)
//    val kernel = KPM.jacksonKernel(order)
//    val rho = range.map(e => KPM.densityOfStates(moments, kernel, e))
//    KPM.plotLines(plot, (range, KPM.integrate(range, rho, moment=1)), "Approx", java.awt.Color.BLACK)
//  }
  
  //scala.tools.nsc.interpreter.ILoop.break(Nil)
  
  var i = 0
  for (f <- dumpdir.listFiles()) {
    implicit val formats = json.DefaultFormats
    val snap = json.Serialization.read[KondoSnap](f.slurp)
    println("t=%g, action=%g".format(snap.time, snap.action))
    drawSpins(snap.spin)
    drawGrid(snap.spin)
    drawGridFft(snap.spin)
//    drawDensity(snap.moments)
    
//    Thread.sleep(500)
    javax.imageio.ImageIO.write(viz.scene.captureImage(), "PNG", new java.io.File("imgs/%03d.png".format(i)))
    //javax.imageio.ImageIO.write(grid.getImage(), "PNG", new java.io.File("imgs2/%03d.png".format(i)))    
    i += 1
  }
}

