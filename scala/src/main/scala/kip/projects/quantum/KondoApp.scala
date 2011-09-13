package kip.projects.quantum


import kip.util.Util.{time, notime}
import ScalarTyp._
import ctor._

object KondoApp extends App {
  val linearSize = 8
  val q = new Quantum(w=linearSize, h=linearSize, t=1, J_eff=2, e_min= -10, e_max= 10)  // hopping only: e_min= -6-0.5, e_max= 3+0.5
  val kpm = new KPM(q.matrix, order=50, nrand=1)
  val e_cutoff = q.scaleEnergy(-2.56)
  val dt = 0.1
  val de = 1e-4
  val fn: R => R = e => if (e < e_cutoff) (e - e_cutoff) else 0
  println("N=%d matrix, %d moments".format(q.matrix.numRows, kpm.order))
  val c = time("Building coefficients. de=%g".format(de))(kpm.expansionCoefficients(de, fn))
  
//  val eig = KPM.eigenvaluesExact(q.matrix)
//  println(eig.size)
//  val idx_cut = eig.size * 3 / 4
//  val eng_cut = (eig(idx_cut-1) + eig(idx_cut)) / 2.0
//  println("Gap between: [%g, %g]".format(eig(idx_cut-1), eig(idx_cut)))
//  println("Choosing potential mu=%g".format(eng_cut))
//  var acc = 0.0
//  for (i <- 0 until idx_cut) {
//    acc += eig(i).re - eng_cut
//  }
//  println(acc)
  
  import kip.graphics._
  import kip.math.Vec3
  val bds = Bounds3d(Vec3(0, 0, 0), Vec3(linearSize-1, linearSize-1, 0))
  val viz = new RetainedScene(bds)
  def drawSpins() {
     val arrows = for (y <- 0 until q.h;
                      x <- 0 until q.w) yield {
      val sx = 0.5*q.field(q.fieldIndex(d=0, x, y))
      val sy = 0.5*q.field(q.fieldIndex(d=1, x, y))
      val sz = 0.5*q.field(q.fieldIndex(d=2, x, y))
      
      val origin = Vec3(x, y, 0)
      val delta  = Vec3(sx.re, sy.re, sz.re)
      new RetainedScene.Arrow(origin, delta, width=0.1)
    }
    
    viz.drawables = Vector(new RetainedScene.Cuboid(bds))
    viz.drawables ++= arrows
    viz.display()
  }
  
//  val cplot = KPM.mkPlot("Coefficients")
//  KPM.plotLines(cplot, (c.indices.toArray.map(i => (i+0.5)/c.size), c.toArray.map(math.abs(_))))
  
  val plot = KPM.mkPlot("Integrated rho")
  for (iter <- 0 until 1000) {
    for (iter2 <- 0 until 20) {
      val r = kpm.randomVector()
      val f0  = notime("Stochastic estimate")(kpm.functionAndGradient(r, c, q.delMatrix))
      q.fieldDerivative(q.delMatrix, q.delField)
      for (i <- q.field.indices) {
        q.field(i) -= dt * q.delField(i)
      }
      q.normalizeField(q.field)
      q.fillMatrix(q.matrix)
      require(math.sqrt((q.matrix - q.matrix.dag).norm2.abs) < 1e-14, "Found non-hermitian hamiltonian!")
    }
    drawSpins()
    
    val eig = notime("Exact diagonalization")(KPM.eigenvaluesExact(q.matrix))
    println("Exact: " + eig.filter(_ < e_cutoff).map(_ - e_cutoff).sum)
    println("Filling: " + eig.filter(_ < e_cutoff).size.toDouble / eig.size)
    KPM.plotLines(plot, (kpm.range, KPM.integrateDeltas(kpm.range, eig, moment=1)), "Exact", java.awt.Color.RED)
    KPM.plotLines(plot, (kpm.range, KPM.integrate(kpm.range, kpm.eigenvaluesApprox(kpm.jacksonKernel), moment=1)), "Approx", java.awt.Color.BLACK)
  }
}
