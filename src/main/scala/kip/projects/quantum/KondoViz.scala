package kip.projects.quantum

import net.liftweb.json
import java.io.File
import scikit.graphics.dim2.Grid
import kip.math.Vec3
import kip.enrich._
import kip.math.fft.FFTComplex
import kip.math.fft.FFTReal
import scala.collection.mutable.ArrayBuffer
import kip.graphics._


object KondoViz extends App {
  val dir = args(0)
  val conf = KondoApp.readConfig(new File(dir+"/cfg.json"))
  
  val w_orig = conf.w
  val h_orig = conf.h
  
  val (left, right, bottom, top) = (0, w_orig, 0, h_orig)
//  val (left, right, bottom, top) = (114, 134, 178, 198) // 'v3/three_quarter/w200/j02_1947_b1_m500_dt20/dump/0022.json'
  
  val w = right-left
  val h = top-bottom
  
  def remapSpinField(field: Array[R]): Array[R] = {
    val ret = new Array[R](3*w*h)
    for (x <- 0 until w; y <- 0 until h) {
      val i = x+w*y
      val xp = (x + left) % w_orig
      val yp = (y + bottom) % h_orig
      val ip = xp + w_orig*yp
      ret(3*i+0) = field(3*ip+0)
      ret(3*i+1) = field(3*ip+1)
      ret(3*i+2) = field(3*ip+2)
    }
    ret
  }
  
  val dumpdir = new java.io.File(dir+"/dump")
  require(dumpdir.isDirectory(), "Cannot load directory %s".format(dumpdir))
  val imgdir = new java.io.File(dir+"/imgs")
  kip.util.Util.createEmptyDir(imgdir)
  
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
  
  
  val latDel1 = Vec3(0.5, -0.5*math.sqrt(3), 0)
  val latDel2 = Vec3(0.5, +0.5*math.sqrt(3), 0)
  val latDel3 = latDel1 + latDel2

  val lat0 = Vec3(0, 0.5*math.sqrt(3)*(h-1), 0)
  val latN = lat0 + latDel1*(h-1) + latDel3*(w-1)

  val bds = Bounds3d(lat0, latN)
  val viz = new RetainedScene(bds, sizew=1200, sizeh=800, cameraDistance=0.9)
  
  val so3 = new KondoSO3()

  
  def readSpin(x: Int, y: Int, field: Array[R]): Vec3 = {
    require(x < w && y < h)
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
  
  def spinDir(field: Array[R]): Array[Vec3] = {
    val ret = for (y <- (h-1) to 0 by -1;
                   x <- 0 until w) yield {
      readSpin(x, y, field)
    }
    ret.toArray
  }
  
  val spinSubLattice: Array[Int] = {
    val ret = for (y <- (h-1) to 0 by -1;
                   x <- 0 until w) yield {
      (x%2) + 2*(y%2)
    }
    ret.toArray
  }

  def drawSpins(field: Array[R]) {
    val sd = spinDir(field)
    
    val arrows = for (i <- 0 until w*h) yield {
      val pos = spinPos(i) + Vec3(0, 0, 1)
      val spin = sd(i)
      val delta = spin*2
      val width = 0.2
      
      val red = java.awt.Color.RED
      val green = java.awt.Color.GREEN
      val black = java.awt.Color.BLACK
      val gray = new java.awt.Color(0, 0, 0, 50)

      if (spinSubLattice(i) == 0)
        new RetainedScene.Arrow(pos, delta, width, color1=black, color2=green)
      else
        new RetainedScene.Arrow(pos, delta, width*0, color1=gray, color2=gray)
    }
    viz.drawables ++= arrows
  }

  
  def drawPlaquettes(field: Array[R]) {
    val plaqs = for (y <- (h-2) to 0 by -1) yield {
      val pts  = new ArrayBuffer[Vec3]
      val vals = new ArrayBuffer[Double]
      
      val p1 = readPos(0, (y+0)%h)
      val p3 = readPos(0, (y+1)%h)
      pts += p3
      pts += p1
      
      for (x <- 0 until w-1) {
        //
        // s3 - s4
        //  \ /  \
        //  s1 - s2
        //
        val p2 = readPos((x+1)%w, (y+0)%h)
        val p4 = readPos((x+1)%w, (y+1)%h)
        val s1 = readSpin((x+0)%w, (y+0)%h, field)
        val s2 = readSpin((x+1)%w, (y+0)%h, field)
        val s3 = readSpin((x+0)%w, (y+1)%h, field)
        val s4 = readSpin((x+1)%w, (y+1)%h, field)
        
        pts += p4
        pts += p2
        vals += s1 dot (s4 cross s3)
        vals += s1 dot (s2 cross s4)
      }
      
      val cg = ColorGradient.blueRed(-0.9, +0.9, alpha=1.0)
      val colors = vals.map(cg.interpolate(_))
      new RetainedScene.TriangleStrip(pts.toArray, colors.toArray)
    }
    
    viz.drawables ++= plaqs
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
    def fftForward(a: Array[R]): Array[Double] = {
      val fft = new FFTReal(Array(h, w))
      val ap = fft.allocFourierArray()
      fft.forwardTransform(a.map(_.toDouble), ap)
      fft.uncompressFourierArray(ap)
    }
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
    
//    // three-Q vectors
//    val i1 = 2*(w/2)+0
//    val i2 = 2*((h/2)*w)+0
//    val i3 = 2*((h/2)*w+w/2)+0
//    val sq1 = Vec3(fx(i1), fy(i1), fz(i1)) / (w*h)
//    val sq2 = Vec3(fx(i2), fy(i2), fz(i2)) / (w*h)
//    val sq3 = Vec3(fx(i3), fy(i3), fz(i3)) / (w*h)
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
    if (true || i == 22) {
    implicit val formats = json.DefaultFormats
    val snap = json.Serialization.read[KondoSnap](f.slurp)
    println("t=%g, action=%g".format(snap.time, snap.action))

    val spin = remapSpinField(snap.spin)
    
    viz.drawables = Vector() // Vector(new RetainedScene.Cuboid(bds))
//    drawSpins(spin)
    drawPlaquettes(spin)
    viz.display()

//    so3.drawPath(so3.loop((w/4,w/4), w/4-1).toArray, spin)

//    drawGrid(spin)
    drawGridFft(spin)
//    drawDensity(snap.moments)
    
//    Thread.sleep(500)
    javax.imageio.ImageIO.write(viz.scene.captureImage(), "PNG", new java.io.File(imgdir+"/%03d.png".format(i)))
    //javax.imageio.ImageIO.write(grid.getImage(), "PNG", new java.io.File("imgs2/%03d.png".format(i)))    
    }
    i += 1
  }
}

