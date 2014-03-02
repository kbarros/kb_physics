package kip.projects.quantum.qmd

import java.io.File

import scala.Array.canBuildFrom
import scala.util.Random

import Units._
import kip.enrich.enrichFile
import kip.math.Vec3
import kip.projects.cuda.JCudaWorld
import kip.projects.quantum.kpm.ComplexKPMCpu
import kip.projects.quantum.kpm.ComplexKPMGpuS
import kip.projects.quantum.kpm.KPMUtil
import kip.util.JacksonWrapper.deserialize
import kip.util.JacksonWrapper.serialize
import kip.util.Util



object TbMD extends App {
  case class Conf(T: Double, gamma: Double, dt: Double, dumpEvery: Int, 
                  M: Int, s: Int, model: Map[String, String])
  case class Snap(time: Double, moments: Array[Double], e_lo: Double, e_hi: Double,
                  energy: Double, x: Array[Vec3], v: Array[Vec3])
  
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
  
  val seed = 0 // System.currentTimeMillis().toInt
  val rand = new Random(seed)
  
  // create output directory for spin configurations
  val dumpdir = new java.io.File(dir+"/dump")
  Util.createEmptyDir(dumpdir)
  
  val conf = loadConf(dir+"/cfg.json")
  
  val pot = GoodwinSi
  val fillingFraction = conf.model("filling").toDouble
  val numAtoms = conf.model("numAtoms").toInt
  val r0 = conf.model("r0").toDouble*angstrom
  val lat = new LinearChain(numAtoms, r0)
  val x = lat.initialPositions
  val v = Array.fill(numAtoms)(Vec3.zero)
  
  println(s"${lat.numAtoms} atoms, ${conf.M} moments")
  
  val kpm = try {
    val cworld = new JCudaWorld(deviceIndex=0)
    new ComplexKPMGpuS(cworld)
  } catch {
    case _: Throwable => {
      println("CUDA device not available! Defaulting to CPU implementation.")
      ComplexKPMCpu
    }
  }
  
  var stepCnt = 0
  var dumpCnt = 0
  
  // don't overwrite any existing data
  // require(!(new File(dumpFilename(dumpCnt))).exists, "Refuse to overwrite dump file %s".format(dumpFilename(dumpCnt)))
  
  def calcMomentsAndDump() {
    return
    val tbh = new TbHamiltonian(pot, lat, x)
    val r = KPMUtil.allVectors(tbh.n)
    val fd = kpm.forward(conf.M, r, tbh.H, KPMUtil.energyScale(tbh.H))
    val mu = tbh.findChemicalPotential(fd, fillingFraction)
    val energy  = tbh.energyAtFixedFilling(kpm, fd, mu, fillingFraction, conf.T)
    
    // dump configuration
    val snap = Snap(time=stepCnt*conf.dt/fs,
                    moments=fd.mu,
                    e_lo=fd.es.lo/eV,
                    e_hi=fd.es.hi/eV,
                    energy=energy/eV,
                    x=x, v=v)
    val filename = dumpFilename(dumpCnt)
    println(s"Dumping $filename (t=${stepCnt*conf.dt/fs} fs)")
    kip.util.Util.writeStringToFile(serialize(snap), filename)
    dumpCnt += 1
  }
  
  def timestep(f: Array[Vec3], m: Double, gamma: Double, T: Double, dt: Double, rand: Random) {
    // update momentum at fixed position
    for (i <- 0 until lat.numAtoms) {
      val eta = Vec3(rand.nextGaussian(), rand.nextGaussian(), rand.nextGaussian())
      v(i) += f(i)*(dt/m) // + v(i)*(-dt/gamma) + eta*math.sqrt(2*kB*T*dt/(m*gamma))
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
  
  println(s"# M=${conf.M}, dt=${conf.dt/fs}")
  
  while (true) {
    Util.notime("Langevin dynamics") (for (_ <- 0 until conf.dumpEvery) {
      val tbh = new TbHamiltonian(pot, lat, x)
      val r = KPMUtil.allVectors(tbh.n)
      val fd = kpm.forward(conf.M, r, tbh.H, KPMUtil.energyScale(tbh.H))
      val mu = tbh.findChemicalPotential(fd, fillingFraction)
      val e_pot  = tbh.energyAtFixedFilling(kpm, fd, mu, fillingFraction, conf.T)
      val e_kin = 0.5 * massSi * v.map(_.norm2).sum
      val e_tot = e_pot + e_kin
      println(s"$stepCnt ${x(0).x}, e_tot=$e_tot e_pot=$e_pot e_kin=$e_kin")
      val f = tbh.force(kpm, fd, mu, conf.T)
      timestep(f, massSi, conf.gamma, conf.T, conf.dt, rand)
    })
    calcMomentsAndDump()
  }
}

