package kip.projects.quantum.kpm

import smatrix._
import scala.util.Random
import kip.projects.cuda.JCudaWorld

object Test extends App {
//  testExpansionCoeffs()
//  testKPM1()
//  testKPM1Cuda()
  testKPM2()
  
  def testExpansionCoeffs() {
    val M = 10
    for (Mq <- Seq(1*M, 10*M, 100*M, 1000*M)) {
      val es = new EnergyScale(lo= -1, hi= +1)
      val f = (e: Double) => if (e < 0.12) e else 0
      val c = KPMUtil.expansionCoefficients(M, Mq, f=f, es=es)
      println(s"quad. points=$Mq, c=<${c.take(10).mkString(",")}>")
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
    val Mq = 4*M
    val es: EnergyScale = KPMUtil.energyScale(H)
    val r = KPMUtil.allVectors(n)
    val fd = ComplexKPMCpu.forward(M, Mq, r, H, es)
    val E = ComplexKPMCpu.function(fd, f)
    val dE_dH = ComplexKPMCpu.gradient(fd, f)
    println(s"KPM energy=$E, expected 50")
    println(s"KPM de_dH(0,0)=${dE_dH(0,0)}, expected 10")
  }

  def testKPM1Cuda() {
    import Constructors.complexDbl._
    
    val n = 4
    val H = {
      val ret = sparse(n, n)
      ret(0,0) = 5.0
      ret(1,1) = -1.0
      ret(2,2) = 0.0
      ret(3,3) = 0.0
      ret.toPacked
    }
    
    def f(x: Double) = x*x
    
    val M = 1000
    val Mq = 4*M
    val es = KPMUtil.energyScale(H)
    val r = KPMUtil.allVectors(n)
    val fd = ComplexKPMCpu.forward(M, Mq, r, H, es)
    val E = ComplexKPMCpu.function(fd, f)
    val dE_dH_cpu = ComplexKPMCpu.gradient(fd, f)
    println(s"CPU E=$E dE_dH(0,0)=${dE_dH_cpu(0,0)}")
    
    val cworld = new JCudaWorld(deviceIndex=0)
    val kpm = new ComplexKPMGpuS(cworld)
    var iter = 0
    while (true) {
      val ms = System.currentTimeMillis()
      val fd = kpm.forward(M, Mq, r, H, es)
      val E = kpm.function(fd, f)
      val dE_dH = kpm.gradient(fd, f)
      println(s"E=$E, dE_dH(0,0)=${dE_dH(0,0)}")
      println(s"Iter $iter Elapsed ${(System.currentTimeMillis() - ms)/1000.0}")
      iter += 1
    }
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
    
    def exactEnergy(H: Dense[Scalar.ComplexDbl]) = {
      H.eig._1.toArray.map(_.re).map(f _).sum
    }
    
    val M = 1000
    val Mq = 4*M
    val es: EnergyScale = KPMUtil.energyScale(H)
    
    val s = 10
    val r = Array(
        KPMUtil.allVectors(n),
        KPMUtil.uncorrelatedVectors(n, s, rand),
        KPMUtil.correlatedVectors(n, s, _%s, rand))
    val fd = r.map(ComplexKPMCpu.forward(M, Mq, _, H, es))
    val E = fd.map(ComplexKPMCpu.function(_, f))
    val dE_dH = fd.map(ComplexKPMCpu.gradient(_, f))
    
    println(s"Energy exact=${exactEnergy(H.toDense)} kpm=${E(0)} stoch1=${E(1)} stoch2=${E(2)}")
    
    val eps = 1e-6
    val dH = sparse(n, n)
    dH(it, jt) = eps
    dH(jt, it) = eps
    
    val de_dH_exact = (exactEnergy((H + dH).toDense)-exactEnergy((H - dH).toDense)) / (2*eps)
    println(s"Deriv exact=$de_dH_exact kpm=${dE_dH(0)(it,jt)+dE_dH(0)(jt,it)} stoch1=${dE_dH(1)(it,jt)+dE_dH(1)(jt,it)} stoch2=${dE_dH(2)(it,jt)+dE_dH(2)(jt,it)}")
  }
}
