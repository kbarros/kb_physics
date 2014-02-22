package kip.projects.quantum.kpm

import smatrix._
import scala.util.Random

object Test extends App {
//  testExpansionCoeffs()
//  testKPM1()
  testKPM2()
  
  def testExpansionCoeffs() {
    val order = 10
    for (qp <- Seq(1*order, 10*order, 100*order, 1000*order)) {
      val es = new EnergyScale(lo= -1, hi= +1)
      val f = (e: Double) => if (e < 0.12) e else 0
      val c = KPMUtil.expansionCoefficients(order, quadPts=qp, f=f, es=es)
      println(s"quad. points=$qp, c=<${c.take(10).mkString(",")}>")
    }
    
    val c_exact = Seq(-0.31601, 0.4801, -0.185462, -0.000775683, 0.0255405, 0.000683006,
                      -0.00536911, -0.000313801, 0.000758494, 0.0000425209)
    println(s"mathematica         c=<${c_exact.mkString(",")}>")
  }
  
  def testKPM1() {
    import Constructors.complexDbl._
    
    val n = 4
    val H = {
      val ret = sparse(n, n)
      ret(0,0) = 5.0
      ret(1,1) = -5.0
      ret(2,2) = 0.0
      ret(3,3) = 0.0
      ret.toPacked
    }
    
    def f(x: Double) = x*x  
    
    val M = 1000
    val es: EnergyScale = KPMUtil.energyScale(H)
    val c = KPMUtil.expansionCoefficients(M=M, quadPts=4*M, f=(e=>e*e), es=es)
    val r = KPMUtil.allVectors(n)
    val (e, de_dH) = ComplexKPMCpu.functionAndGradient(c, r, H, es)
    
    println(s"KPM energy=$e, expected 50")
    println(s"KPM de_dH(0,0)=${de_dH(0,0)}, expected 10")
  }
  
  def testKPM2() {
    import Constructors.complexDbl._
    
    val rand = new Random(2)
    val n = 20
    val it = 3
    val jt = 4
    val H = {
      val ret = sparse(n, n)
      ret(it, jt) = 1.2 - 0.3*I
      ret(jt, it) = 1.2 + 0.3*I
      for (i <- 0 until n)
        ret(i, i) = 0.0
      for (iter <- 0 until 4*n) {
        val i = rand.nextInt(n)
        val j = rand.nextInt(n)
        val r = rand.nextGaussian() + rand.nextGaussian() * I
        ret(i, j) += r
        ret(j, i) += r.conj
      }
      ret.toPacked
    }
    
    def f(x: Double) = x*x
    
    def exactEnergy(H: Dense[Cd]) = {
      H.eig._1.toArray.map(_.re).map(f _).sum
    }
    
    val M = 1000
    val es: EnergyScale = KPMUtil.energyScale(H)
    val c = KPMUtil.expansionCoefficients(M=M, quadPts=4*M, f=(e=>e*e), es=es)
    val r = KPMUtil.allVectors(n)
    val (e, de_dH) = ComplexKPMCpu.functionAndGradient(c, r, H, es)
    
    val s = 10
    
    val r2 = KPMUtil.uncorrelatedVectors(n, s, rand)
    val (e2, de_dH2) = ComplexKPMCpu.functionAndGradient(c, r2, H, es)
    val r3 = KPMUtil.correlatedVectors(n, s, _%s, rand)
    val (e3, de_dH3) = ComplexKPMCpu.functionAndGradient(c, r3, H, es)
    
    println(s"Energy exact=${exactEnergy(H.toDense)} kpm=$e stoch1=$e2 stoch2=$e3")
    
    val eps = 1e-6
    val dH = sparse(n, n)
    dH(it, jt) = eps
    dH(jt, it) = eps
    
    val de_dH_exact = (exactEnergy((H + dH).toDense)-exactEnergy((H - dH).toDense)) / (2*eps)
    println(s"Deriv exact=$de_dH_exact kpm=${de_dH(it,jt)+de_dH(jt,it)} stoch1=${de_dH2(it,jt)+de_dH2(jt,it)} stoch2=${de_dH3(it,jt)+de_dH3(jt,it)}")
  }
}
