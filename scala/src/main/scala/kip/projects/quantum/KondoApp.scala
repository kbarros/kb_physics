package kip.projects.quantum


import kip.util.Util.{time, notime}
import ScalarTyp._
import ctor._

object KondoApp extends App {
  val q = new Quantum(w=10, h=10, t=1, J_eff=2, e_min= -10, e_max= 10)  // hopping only: e_min= -6-0.5, e_max= 3+0.5
  val kpm = new KPM(q.matrix, order=200, nrand=1)
  val e_cutoff = 0.0
  val de = 1e-5
  val fn: R => R = e => if (e < e_cutoff) e else 0
  println("N=%d matrix, %d moments".format(q.matrix.numRows, kpm.order))
  val c = time("Building coefficients. de=%g".format(de))(kpm.expansionCoefficients(de, fn))
  
  val plot = KPM.mkPlot()
  for (iter <- 0 until 100) {
    val r = kpm.randomVector()
    val f0  = notime("Stochastic estimate")(kpm.functionAndGradient(r, c, q.delMatrix))
    val eig = notime("Exact diagonalization")(kpm.eigenvaluesExact())
//    println("exact F(e=0)   = " +eig.filter(_ < 0.0).sum)
//    println("exact F(e=0.2) = " +eig.filter(_ < 0.2).sum)
    q.fieldDerivative(q.delMatrix, q.delField)

    for (i <- q.field.indices) {
      q.field(i) -= 2.0 * q.delField(i)
    }
    q.normalizeField(q.field)

    q.fillMatrix(q.matrix)

    KPM.plotLines(plot, (kpm.range, KPM.integrateDeltas(kpm.range, eig, moment=1)), "Exact", java.awt.Color.RED)
    KPM.plotLines(plot, (kpm.range, KPM.integrate(kpm.range, kpm.eigenvaluesApprox(kpm.jacksonKernel), moment=1)), "Approx", java.awt.Color.BLACK)
  }
}
