package kip.projects.quantum.qmd

import scala.math.sqrt
import scala.util.Random
import Units._
import kip.math.Vec3
import java.io.File
import kip.util.JacksonWrapper._
import kip.util.Util
import kip.enrich._
import kip.projects.cuda.JCudaWorld
import kip.projects.quantum.kpm.ComplexKPM
import kip.projects.quantum.kpm.ComplexKPMCpu
import kip.projects.quantum.kpm.ComplexKPMGpuS



object TbMD extends App {
  case class Conf(T: Double, mu: Double, gamma: Double, dt: Double, dumpEvery: Int, 
                  M: Int, nrand: Int, model: Map[String, String])
  case class Snap(time: Double, moments: Array[Double], e_lo: Double, e_hi: Double,
                  energy: Double, filling: Double, x: Array[Vec3], v: Array[Vec3])
  
  def dumpFilename(iter: Int) = dumpdir+"/%04d.json".format(iter)
  
  def loadConf(filename: String): Conf = {
    import kip.util.JacksonWrapper._
    val s = (new File(filename)).slurp
    val s2 = s.split("\n") map { _.split("""//""")(0) } mkString("\n")
    val conf = deserialize[Conf](s2)
    conf.copy(
        T = conf.T*kelvin,
        mu = conf.mu*eV,
        gamma = conf.gamma/fs,
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
    val tbh = new TbHamiltonian(pot, lat, x)
    
    val moments = tbh.moments(kpm, conf.M)
    val energy  = Double.NaN
    val filling = Double.NaN
    // println("  Action  = %.7g".format(action))
    // println("  Filling = %g".format(filling))
    
    // dump configuration
    val snap = Snap(time=stepCnt*conf.dt/fs,
                    moments=moments,
                    e_lo=tbh.energyScale.lo/eV,
                    e_hi=tbh.energyScale.hi/eV,
                    energy=energy/eV,
                    filling=filling,
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
      v(i) += f(i)*(dt/m) + v(i)*(-gamma*dt) + eta*math.sqrt(2*kB*T*gamma*dt/m)
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
    Util.time("Langevin dynamics") (for (_ <- 0 until conf.dumpEvery) {
      val tbh = new TbHamiltonian(pot, lat, x)
      val (e, f) = tbh.energyAndForces(conf.mu, conf.T, kpm, conf.M)
      timestep(f, massSi, conf.gamma, conf.T, conf.dt, rand)
    })
    calcMomentsAndDump()
  }
}

