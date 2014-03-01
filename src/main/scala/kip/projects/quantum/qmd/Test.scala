package kip.projects.quantum.qmd

import kip.math.Vec3
import kip.projects.quantum.kpm.ComplexKPMCpu
import scala.util.Random
import smatrix._
import smatrix.Constructors.complexDbl._
import Units._
import kip.projects.quantum.kpm.EnergyScale
import kip.projects.quantum.kpm.KPMUtil


object Test extends App {
//  testDimer()
  testDimerJoel()
//  testForce()
//  testDensity()
  
  def testDimer() {
    val pot = GoodwinSi
    val r = pot.r0
    val lat = new LinearChain(numAtoms=2, r0=r)
    val delta = Vec3(1.1, 2.3, -1.8).normalize * (-r)
    
    def buildHamiltonian(delta: Vec3) = {
      new TbHamiltonian(pot, lat, Array(Vec3.zero, delta))
    }
    val tbh = buildHamiltonian(delta)
    
    def sqr(x: Double) = x*x
    val eig1 = Array(
      pot.Δsp/2 + pot.hppπ,
      pot.Δsp/2 + pot.hppπ,
      pot.Δsp/2 - pot.hppπ,
      pot.Δsp/2 - pot.hppπ,
      (pot.hssσ - pot.hppσ)/2 + math.sqrt(sqr(pot.Δsp/2-(pot.hssσ+pot.hppσ)/2)+sqr(pot.hspσ)),
      (pot.hssσ - pot.hppσ)/2 - math.sqrt(sqr(pot.Δsp/2-(pot.hssσ+pot.hppσ)/2)+sqr(pot.hspσ)),
      -(pot.hssσ - pot.hppσ)/2 + math.sqrt(sqr(pot.Δsp/2+(pot.hssσ+pot.hppσ)/2)+sqr(pot.hspσ)),
      -(pot.hssσ - pot.hppσ)/2 - math.sqrt(sqr(pot.Δsp/2+(pot.hssσ+pot.hppσ)/2)+sqr(pot.hspσ))
    ).sorted
    val eig2 = tbh.H.toDense.eig._1.toArray.map(_.re).sorted
    
    val eigDev = math.sqrt((eig1 zip eig2).map(e => sqr(e._1-e._2)).sum)
    println(s"Eigenvalue deviation: $eigDev")
  }
  
  
  def testDimerJoel() { 
    val r = 4.2*bohr
    val pot = GoodwinSi
    val lat = new LinearChain(numAtoms=2, r0=r)
    
    val pos = Array(Vec3.zero, Vec3(r, 0, 0))
    val tbh = new TbHamiltonian(pot, lat, pos)
    
    val natoms = lat.numAtoms
    println(s"Pair energy / atom = ${pot.phi(r) / (natoms*rydberg)} (? 0.17715058505537967 rydberg)")
    
    val eig = tbh.H.toDense.eig._1.toArray.map(_.re).sorted
    
    val nspin = 2
    val elecEnergy = nspin*eig.take(natoms*pot.numFilledOrbitalsPerSite).sum
    println(s"Elec. energy / atom = ${elecEnergy / (natoms*rydberg)} (? -0.37133489682938575 rydberg)")
    
    val mu = 3.15*eV
    val T = 0.001*eV / kB
    val M = 2000
    val kpm = ComplexKPMCpu
    val fd = kpm.forward(M, KPMUtil.allVectors(tbh.n), tbh.H, KPMUtil.energyScale(tbh.H))
    val (e, _) = tbh.energyAndForce(kpm, fd, mu, T)
    val e_kpm = e + mu*natoms*pot.numFilledOrbitalsPerSite*nspin
    val n_kpm = tbh.fillingFraction(kpm, fd, mu, T)
    println(s"Kpm : e=${e_kpm/(natoms*rydberg)} (? -0.19418 rydberg) at n=$n_kpm")
  }
  
  def testForce() {
    val r = 4.2*bohr
    val pot = GoodwinSi
    val lat = new LinearChain(numAtoms=2, r0=r)
    val mu = 2.8525*eV
    val T = 0.1 * eV / kB
    val M = 1000
    val kpm = ComplexKPMCpu

    def calc(pos: Vec3) = {
      val x = Array(Vec3.zero, pos)
      val v = Array(Vec3.zero, Vec3.zero)
      val tbh = new TbHamiltonian(pot, lat, x)
      val fd = kpm.forward(M, KPMUtil.allVectors(tbh.n), tbh.H, KPMUtil.energyScale(tbh.H))
      // println(s"Filling ${tbh.fillingFraction(kpm, fd, mu, T)}")
      tbh.energyAndForce(kpm, fd, mu, T)
    }
    
    val pos = Vec3(1.3, -0.2, 0.7).normalize * r
    println(s"force ${calc(pos)._2(1)} (<1.8277, ...>)")
    
    val eps = 1e-4
    def deriv(dir: Vec3): Double = {
      val ep = calc(pos + dir*eps)._1
      val em = calc(pos - dir*eps)._1
      (ep - em) / (2*eps)
    }
    
    val forceDiscrete = -Vec3(deriv(Vec3(1,0,0)), deriv(Vec3(0,1,0)), deriv(Vec3(0,0,1)))
    println(s"fdisc $forceDiscrete (<1.8277, ...>")
  }
  
  def testDensity() {
    val n = 4
    val H = {
      val ret = sparse(n, n)
      ret(0, 0) = 0.5
      ret(1, 1) = -0.5
      ret(2, 2) = 0.0
      ret(3, 3) = 0.0
      ret.toPacked
    }
    val es = new EnergyScale(-1, 1)
//    val es = KPMUtil.energyScale(H)
    val r = KPMUtil.allVectors(n)
    val M = 100
    val quadPts = 4*M
    val fd = ComplexKPMCpu.forward(M, r, H, es)
    val (x, rho) = KPMUtil.densityFunction(fd.gamma, es)
    scikit.util.Commands.plot(x, rho)
    val (xp, irho) = KPMUtil.integratedDensityFunction(fd.gamma, es)
    scikit.util.Commands.plot(xp, irho)

//    println("total density "+KPMUtil.densityProduct(gamma, x=>x, es))
  }
}
