package kip.projects.quantum

import com.codahale.jerkson.Json._
import kip.util.Util
import ctor._
import java.io.File
import scikit.graphics.dim2.Grid
import kip.math.Vec3


case class KondoConf(w: Int, h: Int, t: Double, J_eff: Double, mu: Double,
                     order: Int, de: Double, dt_per_rand: Double, nrand: Int, dumpPeriod: Int,
                     initConf: String)
case class KondoSnap(time: Double, action: Double, filling: Double, eig: Array[Double], spin: Array[Float])


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
  
  var i = 0
  for (f <- dumpdir.listFiles()) {
    val snap = parse[KondoSnap](f)
    println(snap.time + " "+snap.action)
    drawSpins(snap.spin)
    drawGrid(snap.spin)
//    Thread.sleep(200)
//    javax.imageio.ImageIO.write(viz.scene.captureImage(), "PNG", new java.io.File("imgs/%03d.png".format(i)))
//    javax.imageio.ImageIO.write(grid.getImage(), "PNG", new java.io.File("imgs2/%03d.png".format(i)))    
    i += 1
  }
  
  
    
//  val cplot = KPM.mkPlot("Coefficients")
//  KPM.plotLines(cplot, (c.indices.toArray.map(i => (i+0.5)/c.size), c.toArray.map(math.abs(_))))
//  val plot = KPM.mkPlot("Integrated rho")
//  KPM.plotLines(plot, (kpm.range, KPM.integrateDeltas(kpm.range, eig, moment=1)), "Exact", java.awt.Color.RED)
//  KPM.plotLines(plot, (kpm.range, KPM.integrate(kpm.range, kpm.eigenvaluesApprox(kpm.jacksonKernel), moment=1)), "Approx", java.awt.Color.BLACK)

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
  
  val q = new Quantum(w=w, h=h, t=t, J_eff=J_eff, e_min= -10, e_max= 10)
  val kpm = new KPM(q.matrix, order=order, nrand=nrand)
  initConf match {
    case "random" => q.setFieldRandom(q.field, kpm.rand)
    case "allout" => q.setFieldAllOut(q.field)
    case "threeout" => q.setFieldThreeOut(q.field)
  }
  q.fillMatrix(q.matrix)
  val dt = dt_per_rand * nrand
  val fn: R => R = e => if (e < mu) (e - mu) else 0
  println("N=%d matrix, %d moments".format(q.matrix.numRows, kpm.order))
  val c = Util.time("Building coefficients. de=%g".format(de))(kpm.expansionCoefficients(de, fn))
  
  import kip.projects.cuda._
  val cworld = new JCudaWorld(deviceIndex)
  val ckpm = new CuKPM(cworld, q.matrix, order, nrand)
  
  for (iter <- 0 until 1000) {
    Util.time("Iteration "+iter) (for (iter2 <- 0 until dumpPeriod) {
      val r = kpm.randomVector()
//      val f0  = kpm.functionAndGradient(r, c, q.delMatrix)
      val f0  = ckpm.functionAndGradient(r, c, q.delMatrix)
      q.fieldDerivative(q.delMatrix, q.delField)
      for (i <- q.field.indices) {
        q.field(i) -= dt * q.delField(i)
      }
      q.normalizeField(q.field, validate=true)
      q.fillMatrix(q.matrix)
      require(math.sqrt((q.matrix - q.matrix.dag).norm2.abs) < 1e-14, "Found non-hermitian hamiltonian!")
    })
    
//    val eig = Util.time("Exact diagonalization")(KPM.eigenvaluesExact(q.matrix))
//    val action = eig.filter(_ < mu).map(_ - mu).sum
//    val filling = eig.filter(_ < mu).size.toDouble / eig.size
//    println("Action: " + action)
//    println("Filling: " + filling)
//    println()
    
    val time = iter*dumpPeriod*dt
//    val snap = KondoSnap(time=time, action=action, filling=filling, eig=eig, spin=q.field)
    val snap = KondoSnap(time=time, action=0, filling=0, eig=Array(0f), spin=q.field)
    generate(snap, new File(dumpdir+"/%04d.json".format(iter)))
  }
}
