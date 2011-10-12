package kip.projects.quantum

import com.codahale.jerkson.Json._
import kip.util.Util
import ctor._
import java.io.File


case class KondoConf(w: Int, h: Int, t: Double, J_eff: Double, mu: Double, order: Int, de: Double, dt_per_rand: Double, nrand: Int, dumpPeriod: Int)
case class KondoSnap(time: Double, action: Double, filling: Double, eig: Array[Float], spin: Array[Float])


object KondoViz extends App {
  val dir = args(0)
  val conf = parse[KondoConf](new File(dir+"/cfg.json"))
  val dumpdir = new java.io.File(dir+"/dump")
  require(dumpdir.isDirectory(), "Cannot load directory %s".format(dumpdir))
  
  import kip.graphics._
  import kip.math.Vec3
  val bds = Bounds3d(Vec3(0, 0, 0), Vec3(conf.w-1, conf.h-1, 0))
  val viz = new RetainedScene(bds)
  def drawSpins(field: Array[R]) {
     val arrows = for (y <- 0 until conf.h;
                       x <- 0 until conf.w) yield {
      val sx = 0.5*field(0 + x*3 + y*3*conf.w)
      val sy = 0.5*field(1 + x*3 + y*3*conf.w)
      val sz = 0.5*field(2 + x*3 + y*3*conf.w)
      
      val origin = Vec3(x, y, 0)
      val delta  = Vec3(sx.re, sy.re, sz.re)
      new RetainedScene.Arrow(origin, delta, width=0.1)
    }
    viz.drawables = Vector(new RetainedScene.Cuboid(bds))
    viz.drawables ++= arrows
    viz.display()
  }

  for (f <- dumpdir.listFiles()) {
    val snap = parse[KondoSnap](f)
    println(snap.time + " "+snap.action)
    drawSpins(snap.spin)
  }
  
  
//  val cplot = KPM.mkPlot("Coefficients")
//  KPM.plotLines(cplot, (c.indices.toArray.map(i => (i+0.5)/c.size), c.toArray.map(math.abs(_))))
//  val plot = KPM.mkPlot("Integrated rho")
//  KPM.plotLines(plot, (kpm.range, KPM.integrateDeltas(kpm.range, eig, moment=1)), "Exact", java.awt.Color.RED)
//  KPM.plotLines(plot, (kpm.range, KPM.integrate(kpm.range, kpm.eigenvaluesApprox(kpm.jacksonKernel), moment=1)), "Approx", java.awt.Color.BLACK)

}


object KondoApp extends App {
  val dir = args(0)
  val conf = parse[KondoConf](new File(dir+"/cfg.json"))
  import conf._
  
  // create output directory for spin configurations
  val dumpdir = new java.io.File(dir+"/dump")
  Util.createEmptyDir(dumpdir)
  
  val q = new Quantum(w=w, h=h, t=t, J_eff=J_eff, e_min= -10, e_max= 10)
  val kpm = new KPM(q.matrix, order=order, nrand=nrand)
  q.setFieldRandom(q.field, kpm.rand)
  q.fillMatrix(q.matrix)
  val dt = dt_per_rand * nrand
  val fn: R => R = e => if (e < mu) (e - mu) else 0
  println("N=%d matrix, %d moments".format(q.matrix.numRows, kpm.order))
  val c = Util.time("Building coefficients. de=%g".format(de))(kpm.expansionCoefficients(de, fn))
  
  for (iter <- 0 until 1000) {
    Util.time("Iteration "+iter) (for (iter2 <- 0 until dumpPeriod) {
      val r = kpm.randomVector()
      val f0  = kpm.functionAndGradient(r, c, q.delMatrix)
      q.fieldDerivative(q.delMatrix, q.delField)
      for (i <- q.field.indices) {
        q.field(i) -= dt * q.delField(i)
      }
      q.normalizeField(q.field, validate=true)
      q.fillMatrix(q.matrix)
      require(math.sqrt((q.matrix - q.matrix.dag).norm2.abs) < 1e-14, "Found non-hermitian hamiltonian!")
    })
    
    val eig = Util.time("Exact diagonalization")(KPM.eigenvaluesExact(q.matrix))
    val action = eig.filter(_ < mu).map(_ - mu).sum
    val filling = eig.filter(_ < mu).size.toDouble / eig.size
    println("Action: " + action)
    println("Filling: " + filling)
    println()
    
    val time = iter*dumpPeriod*dt
    val snap = KondoSnap(time=time, action=action, filling=filling, eig=eig, spin=q.field)
    generate(snap, new File(dumpdir+"/%04d.json".format(iter)))
  }
}
