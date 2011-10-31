package kip.projects.quantum

import com.codahale.jerkson.Json._
import kip.util.Util
import ctor._
import java.io.File
import scikit.graphics.dim2.Grid
import kip.math.Vec3


case class KondoConf(w: Int, h: Int, t: Double, J_H: Double, mu: Double,
                     order: Int, de: Double, dt_per_rand: Double, nrand: Int, dumpPeriod: Int,
                     initConf: String)
case class KondoSnap(time: Double,
                     action: Double, action_err: Double, filling: Double, filling_err: Double,
                     moments: Array[Float], spin: Array[Float])

object KondoViz extends App {
  val dir = args(0)
  val conf = parse[KondoConf](new File(dir+"/cfg.json"))
  val w = conf.w
  val h = conf.h
  val dumpdir = new java.io.File(dir+"/dump")
  require(dumpdir.isDirectory(), "Cannot load directory %s".format(dumpdir))
  
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
      val s = readSpin(x, y, field) * 0.8
      val origin = Vec3(x, y, 0)
      val delta  = Vec3(s.x, s.y, s.z)
      new RetainedScene.Arrow(origin, delta, width=0.1, color1=java.awt.Color.RED, color2=java.awt.Color.BLUE)
    }
    viz.drawables = Vector(new RetainedScene.Cuboid(bds))
    viz.drawables ++= arrows
    viz.display()
  }

  val grid = new Grid("Order parameter")
  grid.setScale(-1, 1)
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
    
//    for (y <- 0 until 3;
//         x <- 0 until 3) {
//      val s = readSpin(x, y, field).normalize
//      println("%d %d %s %s %s".format(x, y, s.x, s.y, s.z))
//    }
//    println()
  }
  
  val plot = KPM.mkPlot("Integrated rho")
  def drawDensity(moments: Array[R]) {
    val order = moments.size
    val range = KPM.range(5*order)
    val kernel = KPM.jacksonKernel(order)
    val rho = range.map(e => KPM.densityOfStates(moments, kernel, e))
    KPM.plotLines(plot, (range, KPM.integrate(range, rho, moment=1)), "Approx", java.awt.Color.BLACK)
  }
  
  var i = 0
  for (f <- dumpdir.listFiles()) {
    val snap = parse[KondoSnap](f)
    println("t=%g, action=%g+-%g".format(snap.time, snap.action, snap.action_err))
    drawSpins(snap.spin)
    drawGrid(snap.spin)
//    drawDensity(snap.moments)
    
//    Thread.sleep(200)
//    javax.imageio.ImageIO.write(viz.scene.captureImage(), "PNG", new java.io.File("imgs/%03d.png".format(i)))
//    javax.imageio.ImageIO.write(grid.getImage(), "PNG", new java.io.File("imgs2/%03d.png".format(i)))    
    i += 1
  }
}


object KondoApp extends App {
  if (args.size != 2) {
    println("KondoApp requires <dir> and <device> parameters")
    sys.exit
  }
  val dir = args(0)
  val deviceIndex = args(1).toInt
  val conf = parse[KondoConf](new File(dir+"/cfg.json"))
  import conf._
  
  // create output directory for spin configurations
  val dumpdir = new java.io.File(dir+"/dump")
  Util.createEmptyDir(dumpdir)
  
  val q = new Quantum(w=w, h=h, t=t, J_H=J_H, e_min= -10, e_max= 10)
  val kpm = new KPM(q.matrix, nrand=nrand)
  initConf match {
    case "random" => q.setFieldRandom(q.field, kpm.rand)
    case "allout" => q.setFieldAllOut(q.field)
    case "threeout" => q.setFieldThreeOut(q.field)
  }
  q.fillMatrix(q.matrix)
  val dt = dt_per_rand * nrand
  val mup = q.scaleEnergy(mu)
  
  val fn_action:  (R => R) = e => if (e < mup) (e - mup) else 0
  val fn_filling: (R => R) = e => if (e < mup) (1.0 / q.matrix.numRows) else 0
  val c = KPM.expansionCoefficients(order, de, fn_action)
  
  println("N=%d matrix, %d moments".format(q.matrix.numRows, order))
  
  val ckpm = try {
    import kip.projects.cuda._
    new CuKPM(new JCudaWorld(deviceIndex), q.matrix, nrand)
  } catch {
    case _ => { println("CUDA not found"); null }
  }
  
  for (iter <- 0 until 1000) {
    Util.time("Iteration "+iter) (for (iter2 <- 0 until dumpPeriod) {
      val r = kpm.randomVector()
      val f0 = {
        if (ckpm != null)
          ckpm.functionAndGradient(r, c, q.delMatrix)
        else
          kpm.functionAndGradient(r, c, q.delMatrix)
      }
      q.fieldDerivative(q.delMatrix, q.delField)
      for (i <- q.field.indices) {
        q.field(i) -= dt * q.delField(i)
      }
      q.normalizeField(q.field, validate=true)
      q.fillMatrix(q.matrix)
      require(math.sqrt((q.matrix - q.matrix.dag).norm2.abs) < 1e-14, "Found non-hermitian hamiltonian!")
    })

//    { val eig = Util.time("Exact diagonalization")(KPM.eigenvaluesExact(q.matrix))
//      val action = eig.filter(_ < mup).map(_ - mup).sum
//      val filling = eig.filter(_ < mup).size.toDouble / eig.size
//      println("exact action=%g filling=%g\n".format(action, filling)) }
    
    val (moments, action, action_err, filling, filling_err) = Util.time("Precision moments") {
      val alpha = 4
      val momentsOrder = alpha*order
      val nsamples = dumpPeriod/alpha

      val c_action  = KPM.expansionCoefficients(momentsOrder, de, fn_action)
      val c_filling = KPM.expansionCoefficients(momentsOrder, de, fn_filling)
      val actions  = new Array[Double](nsamples)
      val fillings = new Array[Double](nsamples)
      val moments = Array.fill[R](momentsOrder)(0)
      
      for (iter2 <- 0 until nsamples) {
        val r = kpm.randomVector()
        val momentsOne = {
          if (ckpm != null)
            ckpm.momentsStochastic(momentsOrder, r)
          else
            kpm.momentsStochastic(momentsOrder, r)._1
        }
        
        actions(iter2)  = (c_action,  momentsOne).zipped.map(_*_).sum
        fillings(iter2) = (c_filling, momentsOne).zipped.map(_*_).sum
        for (i <- moments.indices)
          moments(i) += momentsOne(i) / nsamples
      }
      
      import kip.util.Statistics._
      (moments, mean(actions), mean_err(actions), mean(fillings), mean_err(fillings))
    }
    
    println("Action = %g +- %g".format(action, action_err))
    println("Filling = %g +- %g".format(filling, filling_err))
    println()

    val time = iter*dumpPeriod*dt
//    val snap = KondoSnap(time=time, action=action, filling=filling, eig=eig, spin=q.field)
    val snap = KondoSnap(time=time, action=action, action_err=action_err, filling=filling, filling_err=filling_err, moments=moments, spin=q.field)
    generate(snap, new File(dumpdir+"/%04d.json".format(iter)))
  }
}
