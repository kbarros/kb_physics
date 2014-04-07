package kip.projects.quantum.qmd2

import java.io.File
import scala.util.Random
import kip.projects.quantum.qmd.{Potential, GoodwinSi, Lattice, LinearChain, SquareLattice, DiamondLattice}
import kip.projects.quantum.qmd.Units._
import kip.enrich.enrichFile
import kip.math.Vec3
import kip.math.Math.sqr
import kip.projects.cuda.JCudaWorld
import kip.projects.quantum.kpm2.{KPMUtil, KPMComplexCpu, KPMComplexGpu}
import kip.util.JacksonWrapper.deserialize
import kip.util.JacksonWrapper.serialize
import kip.util.Util


object TbMD extends App {
  case class Conf(T: Double, gamma: Double, dt: Double, randomSeed: Int, dumpEvery: Int,
                  M: Int, Mq: Int, s: Int, model: Map[String, String])
  case class Snap(time: Double, energy: Double, mu: Double,
                  bdsLo: Array[Double], bdsHi: Array[Double], 
                  id: Array[Int], x: Array[Double], v: Array[Double],
                  moments: Array[Double], energyScale: (Double, Double))
  
  def dumpFilename(iter: Int) = dumpdir+"/%04d.json".format(iter)
  
  def loadConf(filename: String): Conf = {
    import kip.util.JacksonWrapper._
    val s = (new File(filename)).slurp
    val s2 = s.split("\n") map { _.split("""//""")(0) } mkString("\n")
    val conf = deserialize[Conf](s2)
    conf.copy(
        T = conf.T*kelvin,
        gamma = conf.gamma*fs,
        dt = conf.dt*fs
    )
  }
  
  if (args.size != 2) {
    println("TbMD requires <dir> and <device> parameters")
    sys.exit
  }
  val dir = args(0)
  val deviceIndex = args(1).toInt
    
  // create output directory for spin configurations
  val dumpdir = new java.io.File(dir+"/dump")
  Util.createEmptyDir(dumpdir)
  
  val conf = loadConf(dir+"/cfg.json")
  
  val rand = new Random(conf.randomSeed)
  val pot = GoodwinSi
  val mass = massSi
  val fillingFraction = conf.model("filling").toDouble
  val numAtoms = conf.model("numAtoms").toInt
  val r0 = conf.model("r0").toDouble*angstrom
  
  val periodic = conf.model("periodic").toBoolean
  val lat = conf.model("lattice") match {
    case "linear" => {
      new LinearChain(numAtoms, r0, periodic)
    }
    case "square" => {
      val lx = math.sqrt(numAtoms).round.toInt
      new SquareLattice(lx, lx, r0, periodic)
    }
    case "diamond" => {
      val lx = math.cbrt(numAtoms/8).round.toInt
      new DiamondLattice(lx, lx, lx, r0, periodic)
    }
  }
  val x = lat.initialPositions
  val v = Array.fill(numAtoms)(Vec3.zero)
  val tbh = new TbHamiltonian(pot, lat, x)
  
  println(s"numAtoms=${lat.numAtoms}, M=${conf.M}, T=${conf.T/kelvin} (K), gamma=${conf.gamma/fs} (fs), dt=${conf.dt/fs} (fs)")
  
  val kpm = try {
    val cworld = new JCudaWorld(deviceIndex=0)
    new KPMComplexGpu(cworld, tbh.H, conf.s, conf.M, conf.Mq)
  } catch {
    case _: Throwable => {
      println("CUDA device not available! Defaulting to CPU implementation.")
      new KPMComplexCpu(tbh.H, conf.s, conf.M, conf.Mq)
    }
  }
  
  var stepCnt = 0
  var dumpCnt = 0
  
  // don't overwrite any existing data
  // require(!(new File(dumpFilename(dumpCnt))).exists, "Refuse to overwrite dump file %s".format(dumpFilename(dumpCnt)))
  
  def effectiveTemperature(): Double = {
    def estimateForce() = {
      // kpm.allVectors()
      // kpm.correlatedVectors(tbh.grouping(_, conf.s), rand)
      kpm.uncorrelatedVectors(rand)
      kpm.forward(KPMUtil.energyScale(tbh.H))
      val mu = tbh.findChemicalPotential(kpm, fillingFraction)
      tbh.force(kpm, mu, conf.T)
    }
    tbh.buildHamiltonian()
    // val f1 = estimateForce()
    val f1 = Array.fill(tbh.nAtoms)(Vec3.zero)
    val f2 = estimateForce()
    val dim = 3
    val numEstimates = 1 // 2 if f1 != 0
    val df2 = (f1, f2).zipped.map((f1, f2) => (f1 - f2).norm2).sum / (numEstimates*dim*tbh.nAtoms)
    (df2*conf.dt*conf.gamma) / (2*mass*kB)
  }
  
  def randomizeVectors() {
    if (conf.s < tbh.H.numRows)
      kpm.correlatedVectors(tbh.grouping(_, conf.s), rand)
    else if (conf.s == tbh.H.numRows)
      kpm.allVectors()
    else
      kpm.uncorrelatedVectors(rand)
  }
  
  def calcMomentsAndDump() {
    // println(s"Teff = ${effectiveTemperature()/kelvin} (K)")
    // sys.exit()
    
    tbh.buildHamiltonian()
    randomizeVectors()
    kpm.forward(KPMUtil.energyScale(tbh.H))
    val mu = tbh.findChemicalPotential(kpm, fillingFraction)
    val e_pot = tbh.energyAtFixedFilling(kpm, mu, fillingFraction, conf.T)
    val e_kin = 0.5 * mass * v.map(_.norm2).sum
    val energy = e_pot + e_kin
    
    // compare with exact diagonalization
    if (false) {
      val eig = tbh.H.toSmatrix.toDense.eig._1.toArray.map(_.re).sorted
      val nspin = 2
      val elecEnergy = nspin*eig.take((tbh.n*fillingFraction).round.toInt).sum
      var E = kpm.eval(tbh.localFermiEnergy(_, conf.T, mu))
      val occupiedStates = fillingFraction*nspin*tbh.n
      val elecEnergyKPM = E + mu*occupiedStates
      println(s"exact=$elecEnergy kpm=$elecEnergyKPM")
    }
    
    // dump configuration
    val snap = Snap(time=stepCnt*conf.dt/fs,
                    energy=energy/eV,
                    bdsLo = Array(lat.boundsLow.x, lat.boundsLow.y, lat.boundsLow.z),
                    bdsHi = Array(lat.boundsHigh.x, lat.boundsHigh.y, lat.boundsHigh.z),
                    id=new Array(tbh.nAtoms),
                    x=x.flatMap(r => Array(r.x, r.y, r.z)),
                    v=v.flatMap(r => Array(r.x, r.y, r.z)),
                    mu=mu,
                    moments=kpm.mu,
                    energyScale=(kpm.es.lo/eV, kpm.es.hi/eV))
    val filename = dumpFilename(dumpCnt)
    println(f"Dumping $filename, t=${stepCnt*conf.dt/fs}%g (fs), E=${energy/eV}%g (eV)")
    kip.util.Util.writeStringToFile(serialize(snap), filename)
    dumpCnt += 1
  }
  
  def timestep(f: Array[Vec3], m: Double, gamma: Double, T: Double, dt: Double, rand: Random) {
    // update momentum at fixed position
    for (i <- 0 until lat.numAtoms) {
      val eta = Vec3(rand.nextGaussian(), rand.nextGaussian(), rand.nextGaussian())
      v(i) += f(i)*(dt/m) + v(i)*(-dt/gamma) + eta*math.sqrt(2*kB*T*dt/(m*gamma))
    }
    // update position at fixed momentum
    for (i <- 0 until lat.numAtoms) {
      x(i) += v(i) * dt
    }
    stepCnt += 1
  }
  
  if (dumpCnt == 0) {
    calcMomentsAndDump()
  }
  
  while (true) {
    for (_ <- 0 until conf.dumpEvery) {
      import Util.time
      time("Build hamiltonian")(tbh.buildHamiltonian())
      time("Random vectors")(randomizeVectors())
      val es = time("Energy Scale")(KPMUtil.energyScale(tbh.H))
      time("Forward calc")(kpm.forward(es))
      val mu = time("Chem pot.")(tbh.findChemicalPotential(kpm, fillingFraction))
      val f = time("Force")(tbh.force(kpm, mu, conf.T))
      time("Timestep")(timestep(f, mass, conf.gamma, conf.T, conf.dt, rand))
    }
    calcMomentsAndDump()
  }
}

