package kip.projects.quantum.kpm2

import scala.util.Random
import kip.projects.cuda.JCudaWorld

object Test extends App {
//  testExpansionCoeffs()
  testKPM1()
//  testKPM1Cuda()
  testKPM2()
//  testDensity()
  
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
    val n = 4
    val H = {
      val ret = new SparseCooComplex(n, n)
      ret.add(0, 0,  5.0, 0.0)
      ret.add(1, 1, -5.0, 0.0)
      ret.add(2, 2,  0.0, 0.0)
      ret.add(3, 3,  0.0, 0.0)
      ret.toCsr
    }
    def f(x: Double) = x*x
    val M = 1000
    val Mq = 4*M
    val s = n
    val kpm = new KPMComplexCpu(H, s, M, Mq)
    kpm.allVectors()
    kpm.forward(KPMUtil.energyScale(H))
    val E = kpm.eval(f)
    val dE_dH = kpm.gradient(f)
    println(s"KPM energy=$E, expected 50")
    println(s"KPM de_dH(0,0)=${dE_dH.get_re(0,0)} + I ${dE_dH.get_im(0,0)}, expected 10")
  }
  
  def testKPM1Cuda() {
    val n = 4
    val H = {
      val ret = new SparseCooComplex(n, n)
      ret.add(0, 0,  5.0, 0.0)
      ret.add(1, 1, -1.0, 0.0)
      ret.add(2, 2,  0.0, 0.0)
      ret.add(3, 3,  0.0, 0.0)
      ret.toCsr
    }
    val es = KPMUtil.energyScale(H)
    
    def f(x: Double) = x*x
    
    val M = 1000
    val Mq = 4*M
    val s = n
    val kpm1 = new KPMComplexCpu(H, s, M, Mq)
    kpm1.allVectors()
    kpm1.forward(es)
    val E1 = kpm1.eval(f)
    val dE_dH1 = kpm1.gradient(f)
    println(s"CPU E=$E1 dE_dH(0,0)=${dE_dH1.get_re(0,0)} + I ${dE_dH1.get_im(0,0)}")
    
    val cworld = new JCudaWorld(deviceIndex=0)
    val kpm2 = new KPMComplexGpu(cworld, H, s, M, Mq)
    var iter = 0
    while (true) {
      val ms = System.currentTimeMillis()
      kpm2.allVectors()
      kpm2.forward(es)
      val E2 = kpm2.eval(f)
      val dE_dH2 = kpm2.gradient(f)
      println(s"E=$E2, dE_dH(0,0)=${dE_dH2.get_re(0,0)} + I ${dE_dH2.get_im(0,0)}")
      println(s"Iter $iter Elapsed ${(System.currentTimeMillis() - ms)/1000.0}")
      iter += 1
    }
  }
  
  def testKPM2() {
    // TODO: remove smatrix once sparse add is implemented
    import smatrix._
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
      new SparseCsrComplex(n, n).fromSmatrix(ret.toPacked)
    }
    
    def f(x: Double) = x*x
    
    def exactEnergy(H: Dense[Scalar.ComplexDbl]) = {
      H.eig._1.toArray.map(_.re).map(f _).sum
    }
    
    val M = 1000
    val Mq = 4*M
    val s = n
    val kpm = new KPMComplexCpu(H, s, M, Mq)
    val es: EnergyScale = KPMUtil.energyScale(H)
    
    kpm.allVectors()
    kpm.forward(es)
    val E0 = kpm.eval(f)
    val dEdH0 = kpm.gradient(f).toSmatrix()
    
    kpm.uncorrelatedVectors(rand)
    kpm.forward(es)
    val E1 = kpm.eval(f)
    val dEdH1 = kpm.gradient(f).toSmatrix()
    
    kpm.correlatedVectors(_%s, rand)
    kpm.forward(es)
    val E2 = kpm.eval(f)
    val dEdH2 = kpm.gradient(f).toSmatrix()
    
    println(s"Energy exact=${exactEnergy(H.toSmatrix.toDense)} kpm=${E0} stoch1=${E1} stoch2=${E2}")
    
    val eps = 1e-6
    val dH = sparse(n, n)
    dH(it, jt) = eps
    dH(jt, it) = eps
    
    val de_dH_exact = (exactEnergy((H.toSmatrix + dH).toDense)-exactEnergy((H.toSmatrix - dH).toDense)) / (2*eps)
    println(s"Deriv exact=$de_dH_exact kpm=${dEdH0(it,jt)+dEdH0(jt,it)} stoch1=${dEdH1(it,jt)+dEdH1(jt,it)} stoch2=${dEdH2(it,jt)+dEdH2(jt,it)}")
  }
  
  
    
  def testDensity() {
    val n = 4
    val H = {
      val ret = new SparseCooComplex(n, n)
      ret.add(0, 0,  0.5, 0.0)
      ret.add(1, 1, -0.5, 0.0)
      ret.add(2, 2,  0.0, 0.0)
      ret.add(3, 3,  0.0, 0.0)
      ret.toCsr
    }
    val es = new EnergyScale(-1, 1)
    
    val M = 100
    val Mq = 4*M
    val s = n
    val kpm = new KPMComplexCpu(H, s, M, Mq)
    kpm.allVectors()
    kpm.forward(KPMUtil.energyScale(H))
    
    val (x, rho) = KPMUtil.densityFunction(kpm.gamma, kpm.es)
    scikit.util.Commands.plot(x, rho)
    val (xp, irho) = KPMUtil.integratedDensityFunction(kpm.gamma, kpm.es)
    scikit.util.Commands.plot(xp, irho)
  }
}
