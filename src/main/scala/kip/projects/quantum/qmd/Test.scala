package kip.projects.quantum.qmd

import Units._
import kip.math.Vec3
import smatrix._
import smatrix.Constructors.complexDbl._
import scala.util.Random
import kip.projects.quantum.kpm.ComplexKPMCpu


object Test extends App {
  // testDimer()
  // testDimerJoel()
  testForce()
  
  def testDimer() {
    val pot = new GoodwinSi(rcut=10)
    val r = pot.r0
    val lat = new LinearChain(numAtoms=2, spacing=r)
    val delta = Vec3(1.1, 2.3, -1.8).normalize * (-r)
    
    def buildHamiltonian(delta: Vec3) = {
      new TbHamiltonian(pot, lat, Array(Vec3.zero, delta))
    }
    val tbh = buildHamiltonian(delta)
    
    def sqr(x: Double) = x*x
    val eig1 = Array(
      pot.Δsp + pot.hppπ,
      pot.Δsp + pot.hppπ,
      pot.Δsp - pot.hppπ,
      pot.Δsp - pot.hppπ,
      (pot.hssσ - pot.hppσ)/2 + pot.Δsp/2 + math.sqrt(sqr(pot.Δsp/2-(pot.hssσ+pot.hppσ)/2)+sqr(pot.hspσ)),
      (pot.hssσ - pot.hppσ)/2 + pot.Δsp/2 - math.sqrt(sqr(pot.Δsp/2-(pot.hssσ+pot.hppσ)/2)+sqr(pot.hspσ)),
      -(pot.hssσ - pot.hppσ)/2 + pot.Δsp/2 + math.sqrt(sqr(pot.Δsp/2+(pot.hssσ+pot.hppσ)/2)+sqr(pot.hspσ)),
      -(pot.hssσ - pot.hppσ)/2 + pot.Δsp/2 - math.sqrt(sqr(pot.Δsp/2+(pot.hssσ+pot.hppσ)/2)+sqr(pot.hspσ))
    ).sorted
    val eig2 = tbh.H.toDense.eig._1.toArray.map(_.re).sorted
    
    val eigDev = math.sqrt((eig1 zip eig2).map(e => sqr(e._1-e._2)).sum)
    println(s"Eigenvalue deviation: $eigDev")
    
    val rand = new Random(0)
    val rmat = tbh.H.duplicate.transform(x => if (x == 0) 0 else rand.nextGaussian())
    val de_dh = (rmat + rmat.dag).toPacked
    
    // fake energy with simple, random dependence on elements of H 
    def energy(H: PackedSparse[Scalar.ComplexDbl]) = (de_dh.toDense dagDot H.toDense).re
    
    val dx = 0.0001
    def deriv(dir: Vec3): Double = {
      val ep = energy(buildHamiltonian(delta + dir*dx).H)
      val em = energy(buildHamiltonian(delta - dir*dx).H)
      (ep - em) / (2*dx)
    }
    val forceDiscrete = -Vec3(deriv(Vec3(1,0,0)), deriv(Vec3(0,1,0)), deriv(Vec3(0,0,1)))
    val forceAnalytical = tbh.forces(de_dh)(1)
    
    // println(s"$forceDiscrete $forceAnalytical")
    println(s"Force deviation: ${forceAnalytical-forceDiscrete}")
  }
  
  
  def testDimerJoel() { 
    val r = 4.2*bohr
    val pot = new GoodwinSi(rcut=10)
    val lat = new LinearChain(numAtoms=2, spacing=r)
    
    val pos = Array(Vec3.zero, Vec3(r, 0, 0))
    val tbh = new TbHamiltonian(pot, lat, pos)
    
    val natoms = 2
    println(s"Pair energy / atom = ${pot.phi(r) / (natoms*rydberg)} (? 0.17715058505537967 rydberg)")
    
    val eig = tbh.H.toDense.eig._1.toArray.map(_.re).sorted
    
    val ezero = natoms*8.295*eV
    val nspin = 2
    val elecEnergy = nspin * (eig.take(4).sum - ezero)
    println(s"Elec. energy / atom = ${elecEnergy / (natoms*rydberg)} (? -0.37133489682938575 rydberg)")
  }
  
  def testForce() {
    val r = 4.2*bohr
    val pot = new GoodwinSi(rcut=10)
    val lat = new LinearChain(numAtoms=2, spacing=r)
    val T = 0.1 * eV / kB
    val M = 10
    val mu = 7.0
    // val mu = 7.287906987689409 // half filling
    
    println("M = "+M)
    
    def calc(pos: Vec3) = {
      val x = Array(Vec3.zero, pos)
      val v = Array(Vec3.zero, Vec3.zero)
      val tbmd = new TbMD(x, v, pot, lat, ComplexKPMCpu, M)
      tbmd.energyAndForces(mu, T)
    }
    
    val pos = Vec3(1.3, -0.2, 0.7).normalize * r
    println("energy " + calc(pos)._1)
    println("force " + calc(pos)._2(1))
    
    val eps = 1e-4
    def deriv(dir: Vec3): Double = {
      val ep = calc(pos + dir*eps)._1
      val em = calc(pos - dir*eps)._1
      (ep - em) / (2*eps)
    }
    
    val forceDiscrete = -Vec3(deriv(Vec3(1,0,0)), deriv(Vec3(0,1,0)), deriv(Vec3(0,0,1)))
    println("fdisc "+forceDiscrete)
  }
}
