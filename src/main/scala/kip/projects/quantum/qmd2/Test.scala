package kip.projects.quantum.qmd2


import kip.math.Vec3
import scala.util.Random
import kip.projects.quantum.qmd.{Potential, GoodwinSi, Lattice, LinearChain}
import kip.projects.quantum.qmd.Units._
import kip.projects.quantum.kpm2.KPMComplex
import kip.projects.quantum.kpm2.KPMUtil
import kip.projects.quantum.kpm2.EnergyScale
import kip.projects.quantum.kpm2.SparseCsrComplex
import kip.projects.quantum.kpm2.SparseCooComplex
import kip.projects.quantum.kpm2.KPMComplexCpu


object Test extends App {
  testDimer()
  testDimerJoel()
  testForce()
  
  def testDimer() {
    val pot = GoodwinSi
    val r = pot.r0
    val lat = new LinearChain(numAtoms=2, r0=r, false)
    val delta = Vec3(1.1, 2.3, -1.8).normalize * (-r)
    
    val tbh = new TbHamiltonian(pot, lat, Array(Vec3.zero, delta))
    tbh.buildHamiltonian()
    
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
    val eig2 = tbh.H.toSmatrix.toDense.eig._1.toArray.map(_.re).sorted
    
    val eigDev = math.sqrt((eig1 zip eig2).map(e => sqr(e._1-e._2)).sum)
    println(s"Eigenvalue deviation: $eigDev")
  }
  
  
  def testDimerJoel() { 
    val r = 4.2*bohr
    val pot = GoodwinSi
    val fillingFraction = 0.5
    val T = 0.0
    val lat = new LinearChain(numAtoms=2, r0=r, false)
    val x = lat.initialPositions()
    val tbh = new TbHamiltonian(pot, lat, x)
    tbh.buildHamiltonian()
    
    val natoms = lat.numAtoms
    println(s"Pair energy / atom = ${pot.phi(r) / (natoms*rydberg)} (? 0.17715058505537967 rydberg)")
    
    val eig = tbh.H.toSmatrix.toDense.eig._1.toArray.map(_.re).sorted
    val nspin = 2
    val elecEnergy = nspin*eig.take((tbh.n*fillingFraction).round.toInt).sum
    println(s"Elec. energy / atom = ${elecEnergy / (natoms*rydberg)} (? -0.37133489682938575 rydberg)")
    
    println(s"Total = ${(pot.phi(r)+elecEnergy) / (natoms*rydberg)} rydberg")
    
    val kpm = new KPMComplexCpu(tbh.H, s=tbh.n, M=500, Mq=2000)
    kpm.allVectors()
    kpm.forward(KPMUtil.energyScale(tbh.H))
    val mu = tbh.findChemicalPotential(kpm, fillingFraction)
    val E_kpm = tbh.energyAtFixedFilling(kpm, mu, fillingFraction, T)
    val n_kpm = tbh.fillingFraction(kpm, mu, T)
    println(s"Kpm : e=${E_kpm/(natoms*rydberg)} (? -0.19418 rydberg) at n=$n_kpm ($fillingFraction)")
  }
  
  def testForce() {
    val r = 4.2*bohr
    val pot = GoodwinSi
    val mu = 2.8525*eV
    val T = 0.1 * eV / kB
    val lat = new LinearChain(numAtoms=2, r0=r, false)
    val x = lat.initialPositions()
    val tbh = new TbHamiltonian(pot, lat, x)
    val kpm = new KPMComplexCpu(tbh.H, s=tbh.n, M=1000, Mq=4000)    
    
    def calc(pos: Vec3) = {
      x(1) = pos
      tbh.buildHamiltonian()
      kpm.allVectors()
      kpm.forward(KPMUtil.energyScale(tbh.H))
      (tbh.energy(kpm, mu, T), tbh.force(kpm, mu, T))
    }
    
    val pos = Vec3(1.3, -0.2, 0.7).normalize * r
    println(s"force ${calc(pos)._2(1)} (? <2.038, ...>)")
    
    val eps = 1e-4
    def deriv(dir: Vec3): Double = {
      val ep = calc(pos + dir*eps)._1
      val em = calc(pos - dir*eps)._1
      (ep - em) / (2*eps)
    }
    
    val forceDiscrete = -Vec3(deriv(Vec3(1,0,0)), deriv(Vec3(0,1,0)), deriv(Vec3(0,0,1)))
    println(s"fdisc $forceDiscrete")
  }
}
