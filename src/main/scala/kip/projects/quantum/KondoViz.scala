package kip.projects.quantum

import java.io.File
import scikit.graphics.dim2.Grid
import kip.math.Vec3
import kip.enrich._
import kip.math.fft.FFTComplex
import kip.math.fft.FFTReal
import scala.collection.mutable.ArrayBuffer
import kip.graphics._

// import kip.projects.quantum.prb13.KondoSO3


object KondoViz extends App {
  val dir = args(0)
  val conf = KondoApp.loadKondoConf(dir+"/cfg.json")
  val q = KondoHamiltonian.fromMap(conf.model)
  
  val dumpdir = new java.io.File(dir+"/dump")
  require(dumpdir.isDirectory(), "Cannot load directory %s".format(dumpdir))
  val imgdir = new java.io.File(dir+"/imgs")
  kip.util.Util.createEmptyDir(imgdir)
  
  val bds = {
	val (xs, ys, zs) = q.latticePositions.map(v => (v.x, v.y, v.z)).unzip3
    Bounds3d(Vec3(xs.min, ys.min, zs.min), Vec3(xs.max, ys.max, zs.max))
  }
  
  val viz = new RetainedScene(bds, sizew=1200, sizeh=800, cameraDistance=0.9)
  
//  def triangularSubLattice(i: Int) = {
//    val w = q.asInstanceOf[TriangularLattice].w
//    val x = i % w
//    val y = i / w
//    (x%2) + 2*(y%2)
//  }
//  
//  def kagomeSubLattice(i: Int) = i % 3
  
  def drawSpins(field: Array[R]) {    
    val arrows = for (i <- 0 until q.numLatticeSites) yield {
      val spin = Vec3(field(3*i+0), field(3*i+1), field(3*i+2))
      val pos = q.latticePositions(i) - spin*0.3
      val delta = spin
      val width = 0.1
      
      val red = java.awt.Color.RED
      val green = java.awt.Color.GREEN
      val black = java.awt.Color.BLACK
      val gray = new java.awt.Color(0, 0, 0, 50)

//      if (spinSubLattice(i) == 0)
//        new RetainedScene.Arrow(pos, delta, width, color1=black, color2=green)
//      else
//        new RetainedScene.Arrow(pos, delta, width*0, color1=gray, color2=gray)
      new RetainedScene.Arrow(pos, delta, width, color1=black, color2=green)
    }
    viz.drawables ++= arrows
  }
  
  
  def drawChiralPlaquettes(field: Array[R]) {
    val cg = ColorGradient.blueRed(-0.9, +0.9, alpha=1.0)
    
    q match {
      case q: TriangularLattice => {
        def pos(x: Int, y: Int): Vec3  = q.latticePositions(q.coord2idx(x%q.w, y%q.h))
        def spin(x: Int, y: Int): Vec3 = {
          val i = q.coord2idx(x%q.w, y%q.h)
          Vec3(field(3*i+0), field(3*i+1), field(3*i+2))
        }
        val plaqs = for (y <- (q.h-2) to 0 by -1) yield {
          val pts  = new ArrayBuffer[Vec3]
          val vals = new ArrayBuffer[Double]
          val p1 = pos(0, y+0)
          val p3 = pos(0, y+1)
      	  pts += p3
      	  pts += p1
    	    for (x <- 0 until q.w-1) {
    	      //
      	    // s3 - s4
      	    //  \ /  \
      	    //  s1 - s2
      	    //
      	    val p2 = pos(x+1, y+0)
      	    val p4 = pos(x+1, y+1)
      	    val s1 = spin(x+0, y+0)
      	    val s2 = spin(x+1, y+0)
      	    val s3 = spin(x+0 ,y+1)
      	    val s4 = spin(x+1, y+1)
      	    pts += p4
      	    pts += p2
      	    vals += s1 dot (s4 cross s3)
      	    vals += s1 dot (s2 cross s4)
      	  }
          val colors = vals.map(cg.interpolate(_))
          new RetainedScene.TriangleStrip(pts.toArray, colors.toArray)
        }
        viz.drawables ++= plaqs
      }
      
      case q: KagomeLattice => {
        def posAndSpin(v: Int, x: Int, y: Int): (Vec3, Vec3) = {
          val i = q.coord2idx(v, (x+q.w)%q.w, (y+q.h)%q.h)
          val p = q.latticePositions(i)
          val s = Vec3(field(3*i+0), field(3*i+1), field(3*i+2))
          (p, s)
        }
        
        for (y <- 0 until q.h;
             x <- 0 until q.w) {
          //     D\       /E
          //     - 0 --- 2 - 
          //         \ /  
          //          1   
          //         /B\  
          //     - 2 --- 0 -
          val (pB0, sB0) = posAndSpin(0, x, y)
          val (pB1, sB1) = posAndSpin(1, x, y)
          val (pB2, sB2) = posAndSpin(2, x, y)
          val (pE2, sE2) = posAndSpin(2, x, y+1)
          val (pD0, sD0) = posAndSpin(0, x-1, y+1)
          
          val chi1 = sB0 dot (sB1 cross sB2)
          val tri1 = new RetainedScene.Triangles(Array(pB0, pB1, pB2), cg.interpolate(chi1))
          viz.drawables :+= tri1
          
          if (y < q.h-1 && x > 0) {
            val chi2 = sE2 dot (sD0 cross sB1)
            val tri2 = new RetainedScene.Triangles(Array(pE2, pD0, pB1), cg.interpolate(chi2))
            viz.drawables :+= tri2
          }
        }
      }
    }
  }
  
//  val grid = new Grid("Order parameter")
//  grid.setScale(-0.77, 0.77)
//  scikit.util.Utilities.frame(grid.getComponent(), grid.getTitle())
//  val gridData = new Array[Double]((2*w)*(2*h))
//  def drawGrid(field: Array[R]) {
//    for (y <- 0 until h;
//         x <- 0 until w) yield {
//      //
//      // s4 - s3
//      //  | / |
//      // s1 - s2
//      //
//      val s1 = readSpin((x+0)%w, (y+0)%h, field)
//      val s2 = readSpin((x+1)%w, (y+0)%h, field)
//      val s3 = readSpin((x+1)%w, (y+1)%h, field)
//      val s4 = readSpin((x+0)%w, (y+1)%h, field)
//      gridData((2*y+0)*(2*w) + 2*x+0) = s1 dot (s2 cross s3)
//      gridData((2*y+1)*(2*w) + 2*x+0) = s1 dot (s3 cross s4)
//      gridData((2*y+0)*(2*w) + 2*x+1) = s1 dot (s2 cross s3)
//      gridData((2*y+1)*(2*w) + 2*x+1) = s1 dot (s3 cross s4)
//    }
//    grid.registerData(2*w, 2*h, gridData)
//  }

//  // Fourier transform
//  def splitFieldComponents(field: Array[R]): Seq[Array[R]] = {
//    val x = new Array[R](w*h)
//    val y = new Array[R](w*h)
//    val z = new Array[R](w*h)
//    for (i <- 0 until w*h) {
//      x(i) = field(3*i+0)
//      y(i) = field(3*i+1)
//      z(i) = field(3*i+2)
//    }
//    Seq(x, y, z)
//  }
//  val gridFft = new Grid("FFT")
//  scikit.util.Utilities.frame(gridFft.getComponent(), gridFft.getTitle())
//  val gridFftData = new Array[Double](w*h)
//  def drawGridFft(field: Array[R]) {
//    def fftForward(a: Array[R]): Array[Double] = {
//      val fft = new FFTReal(Array(h, w))
//      val ap = fft.allocFourierArray()
//      fft.forwardTransform(a.map(_.toDouble), ap)
//      fft.uncompressFourierArray(ap)
//    }
//    val Seq(sx, sy, sz) = splitFieldComponents(field)
//    val fx = fftForward(sx)
//    val fy = fftForward(sy)
//    val fz = fftForward(sz)
//    gridFftData.transform(_ => 0)
//    for (i <- 0 until h*w) {
//      import kip.math.Math.sqr
//      for (f <- Seq(fx, fy, fz)) {
//        gridFftData(i) += sqr(f(2*i+0)) + sqr(f(2*i+1))
//      }
//    }
//    gridFft.registerData(w, h, gridFftData)
//    
////    // three-Q vectors
////    val i1 = 2*(w/2)+0
////    val i2 = 2*((h/2)*w)+0
////    val i3 = 2*((h/2)*w+w/2)+0
////    val sq1 = Vec3(fx(i1), fy(i1), fz(i1)) / (w*h)
////    val sq2 = Vec3(fx(i2), fy(i2), fz(i2)) / (w*h)
////    val sq3 = Vec3(fx(i3), fy(i3), fz(i3)) / (w*h)
//  }

  
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
    val snap = kip.util.JacksonWrapper.deserialize[KondoSnap](f.slurp)
    println("t=%g, action=%g".format(snap.time, snap.action))
    
    viz.drawables = Vector() // Vector(new RetainedScene.Cuboid(bds))
    drawSpins(snap.spin)
    drawChiralPlaquettes(snap.spin)
    viz.display()

//    so3.drawPath(so3.loop((w/4,w/4), w/4-1).toArray, spin)

//    drawGrid(spin)
//    drawGridFft(spin)
//    drawDensity(snap.moments)
    
//    Thread.sleep(500)
//    javax.imageio.ImageIO.write(viz.scene.captureImage(), "PNG", new java.io.File(imgdir+"/%03d.png".format(i)))
    //javax.imageio.ImageIO.write(grid.getImage(), "PNG", new java.io.File("imgs2/%03d.png".format(i)))    
    }
    i += 1
  }
}

