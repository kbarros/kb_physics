package kip.projects.quantum

import java.io.File
import scikit.graphics.dim2.Grid
import kip.math.Vec3
import kip.enrich._
import kip.math.fft.FFTComplex
import kip.math.fft.FFTReal
import scala.collection.mutable.ArrayBuffer
import kip.graphics._
import java.awt.event.MouseEvent
import kip.math.Quaternion

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
  
  // independently rotation spins
  var spinRotation = Quaternion.fromAxisAngle(0, 0, 0)
  val mouseHandler = new DragHandler {
    def mouseDraggedDelta(dx: Int, dy: Int, e: MouseEvent) {
      if (e.isAltDown()) {
        val dpp = 0.003 // dimensionless displacement per pixel
        val rpp = 0.01  // radians per pixel
        val q = Quaternion.fromAxisAngle(dy*rpp, dx*rpp, 0)
        spinRotation = (q * spinRotation).normalize
        drawScene()
      }
    }
  }
  viz.scene.canvas.addMouseListener(mouseHandler)
  viz.scene.canvas.addMouseMotionListener(mouseHandler)
  
//  def triangularSubLattice(i: Int) = {
//    val w = q.asInstanceOf[TriangularLattice].w
//    val x = i % w
//    val y = i / w
//    (x%2) + 2*(y%2)
//  }
//  
  def kagomeSubLattice(i: Int) = {
    val q2 = q.asInstanceOf[KagomeLattice]
    val (v, x, y) = q2.idx2coord(i)
    v // v + ((x%2) + 2*(y%2))
  }
  
  def drawSpins(field: Array[R]) {    
    val arrows = for (i <- 0 until q.numLatticeSites) yield {
      val spin = spinRotation.rotate(Vec3(field(3*i+0), field(3*i+1), field(3*i+2)))
      val pos = q.latticePositions(i) - spin*0.3
      val delta = spin*0.85
      val width = 0.1
      
      val red = java.awt.Color.RED
      val green = java.awt.Color.GREEN
      val black = java.awt.Color.BLACK
      val gray = new java.awt.Color(0, 0, 0, 50)
      
//      if (kagomeSubLattice(i) == 2)
        viz.drawables :+= new RetainedScene.Arrow(pos, delta, width, color1=black, color2=green)
    }
//    viz.drawables ++= arrows
  }
  
  
  def drawChiralPlaquettes(field: Array[R]) {
    val cg = ColorGradient.blueRed(-1.0, +1.0, alpha=0.9)
    
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
        val blue = Array.fill(3)(new java.awt.Color(0.8f, 0.8f, 1.0f))
        
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
          val ls1 = new RetainedScene.LineStrip(Array(pB0, pB1, pB2, pB0), blue)
          // viz.drawables :+= ls1
          viz.drawables :+= tri1          
          
          if (y < q.h-1 && x > 0) {
            val chi2 = sE2 dot (sD0 cross sB1)
            val tri2 = new RetainedScene.Triangles(Array(pE2, pD0, pB1), cg.interpolate(chi2))
            val ls2 = new RetainedScene.LineStrip(Array(pE2, pD0, pB1, pE2), blue)
            // viz.drawables :+= ls2
            viz.drawables :+= tri2
          }
        }
      }
      
      case q: SquareLattice => {
        if (false) {
        for (y <- 0 until q.h;
             x <- 0 until q.w) {
          val i = q.coord2idx(x, y)
          // val Sz = field(3*i+1)
          val Sz = spinRotation.rotate(Vec3(field(3*i+0), field(3*i+1), field(3*i+2))).z
          val a = 1.0
          val p1 = q.latticePositions(i) + Vec3(-0.5, -0.5, 0)*a
          val p2 = p1 + Vec3(a, 0, 0)
          val p3 = p1 + Vec3(0, a, 0)
          val p4 = p1 + Vec3(a, a, 0)
          val c = cg.interpolate(Sz)
          val tris =  new RetainedScene.Triangles(Array(p1, p2, p3, p3, p2, p4), c)
          viz.drawables :+= tris
        }
        }
        else {
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
      }
    }
  }

  // Fourier transform
  lazy val gridFft = {
    val ret = new Grid("FFT")
    ret.setAutoScale()
    scikit.util.Utilities.frame(ret.getComponent(), "FFT")
    ret
  }
  
  def drawGridFft(field: Array[Float]) {
    import kip.math.Math.sqr
    q match {
      case q: SquareLattice => {
        val n = q.numLatticeSites
        def splitFieldComponents(field: Array[Float]): Seq[Array[Float]] = {
          val x = new Array[Float](n)
          val y = new Array[Float](n)
          val z = new Array[Float](n)
          for (i <- 0 until n) {
            x(i) = field(3*i+0)
            y(i) = field(3*i+1)
            z(i) = field(3*i+2)
          }
          Seq(x, y, z)
        }
        val fft = new FFTReal(Array(q.h, q.w))
        def fftForward(a: Array[Float]): Array[Double] = {
          val ap = fft.allocFourierArray()
          fft.forwardTransform(a.map(_.toDouble), ap)
          fft.uncompressFourierArray(ap)
        }
        val Seq(sx, sy, sz) = splitFieldComponents(field)
        val fx = fftForward(sx)
        val fy = fftForward(sy)
        val fz = fftForward(sz)
        val gridFftData = Array.fill(n)(0.0)
        for (i <- 0 until n) {
          for (f <- Seq(fx, fy, fz)) {
            gridFftData(i) += sqr(f(2*i+0)) + sqr(f(2*i+1))
          }
        }
        
        val imax = gridFftData.zipWithIndex.maxBy(_._1)._2
        val kx = math.min(imax%q.w, q.w - imax%q.w).toDouble
        val ky = math.min(imax/q.w, q.h - imax/q.w).toDouble
        val lambda_x = q.w / kx
        val lambda_y = q.h / ky
        val lambda = 1 / math.sqrt(1.0/sqr(lambda_x) + 1.0/sqr(lambda_y))
        println(s"kx=$kx ky=$ky lambda=$lambda diagonal=${lambda/math.sqrt(2.0)}")
        
        gridFft.registerData(q.w, q.h, gridFftData)
      }
      
      case _ => ()
    }
  }

  
  lazy val plot = KPM.mkPlot("Integrated rho")
  def drawDensity(moments: Array[R]) {
    val order = moments.size
    val range = KPM.range(5*order)
    val kernel = KPM.jacksonKernel(order)
    val rho = range.map(e => KPM.densityOfStates(moments, kernel, e) / (2*q.numLatticeSites))
    KPM.plotLines(plot, (range, KPM.integrate(range, rho, moment=0)), "Approx", java.awt.Color.BLACK)
  }
  
  //scala.tools.nsc.interpreter.ILoop.break(Nil)
  
  var snap: KondoSnap = null
  
  def drawScene() {
    if (snap != null) {
      viz.drawables = Vector() // Vector(new RetainedScene.Cuboid(bds))
      val field = if (true) {
        snap.spin
      } else {
        q.asInstanceOf[KagomeLattice].setFieldCoplanar3(q.field)
        q.field
      }
      drawSpins(field)
      drawChiralPlaquettes(field)
      viz.display()
    }
  }
  
  var i = 0
  for (f <- dumpdir.listFiles()) {
    if (true || i == 22) {
    snap = kip.util.JacksonWrapper.deserialize[KondoSnap](f.slurp)
    println(s"t=${snap.time}, action=${snap.action} filling=${snap.filling}")
    
    drawScene()
    
    drawGridFft(snap.spin)
    drawDensity(snap.moments)
    
//    Thread.sleep(500)
//    javax.imageio.ImageIO.write(viz.scene.captureImage(), "PNG", new java.io.File(imgdir+"/%03d.png".format(i)))
    //javax.imageio.ImageIO.write(grid.getImage(), "PNG", new java.io.File("imgs2/%03d.png".format(i)))    
    }
    i += 1
  }
}

