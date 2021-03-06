package kip.projects.quantum

import java.io.File
import kip.util.Util
import kip.util.JacksonWrapper._
import kip.enrich._
import ctor._
import scala.math._


case class KondoConf(T: Double, mu: Double,
		             order: Int, order_exact: Int, dt_per_rand: Double,
                     nrand: Int, dumpPeriod: Int, initSpin: String, model: Map[String, String])

case class KondoSnap(time: Double, action: Double, filling: Double, eig: Array[Double],
                     moments: Array[Float], spin: Array[Float])


object KondoApp extends App {
  
  def dumpFilename(iter: Int) = dumpdir+"/%04d.json".format(iter)

  def loadKondoConf(filename: String): KondoConf = {
    import kip.util.JacksonWrapper._
    val s = (new File(filename)).slurp
    val s2 = s.split("\n") map { _.split("""//""")(0) } mkString("\n")
    deserialize[KondoConf](s2)
  }

  if (args.size != 2) {
    println("KondoApp requires <dir> and <device> parameters")
    sys.exit
  }
  val dir = args(0)
  val deviceIndex = args(1).toInt
  
  val seed = 0 // System.currentTimeMillis().toInt
  
  // create output directory for spin configurations
  val dumpdir = new java.io.File(dir+"/dump")
  Util.createEmptyDir(dumpdir)
  
  val conf = loadKondoConf(dir+"/cfg.json")
  import conf._
  
  val q = KondoHamiltonian.fromMap(model)
  
  val kpm = try {
    import kip.projects.cuda._
    new CuKPM(new JCudaWorld(deviceIndex), q.matrix, nrand, seed)
  } catch {
    case _: Throwable => {
      println("CUDA device not available! Defaulting to CPU implementation.")
      new KPM(q.matrix, nrand, seed)
    }
  }
   
  var iter = 0
  
  initSpin match {
    case "random" => {
      q.setFieldRandom(q.field, kpm.rand)
    }
    case s if s.startsWith("dump") => {
      val lastDump = s.replace("dump", "").toInt
      val f = new File(dumpFilename(lastDump))
      println("Reading initial configuration from file: "+f)
      val snap = deserialize[KondoSnap](f.slurp)
      snap.spin.copyToArray(q.field)
      iter = lastDump+1
    }
    case s => {
      q.setField(s)
    }
  }

  // don't overwrite any existing data
  require(!(new File(dumpFilename(iter))).exists, "Refuse to overwrite dump file %s".format(dumpFilename(iter)))
  
  q.fillMatrix(q.matrix)
  val dt = dt_per_rand * nrand
  val mup = q.scaleEnergy(mu)
  val Tp = q.scaleEnergy(T)
  
  // calculates the action in **scaled coordinates**
  // this means that action and its derivative will be too small by a factor of q.e_scale (typically 10)
  val fn_action:  (R => R) = {
    if (Tp == 0) (e => if (e - mup < 0) (e - mup) else 0)
    else { e =>
      val x: Double = (e - mup)/Tp
      if (x < -20)
        (e - mup)
      else if (x > 20)
        0
      else
        -Tp*log(1 + exp(-x))
    }
  }
  
  val fn_filling: (R => R) = e => if (e < mup) (1.0 / q.matrix.numRows) else 0
  
  val c_action   = KPM.expansionCoefficients2(order, quadPts=10*order, fn_action)
  val c2_action  = KPM.expansionCoefficients2(order_exact, quadPts=10*order, fn_action)
  val c2_filling = KPM.expansionCoefficients2(order_exact, quadPts=10*order, fn_filling)
  
  println("N=%d matrix, %d moments".format(q.matrix.numRows, order))
  
  
  def calcMomentsAndDump() {
    // exact moments
    q.fillMatrix(q.matrix)
    val (eig, moments, action, filling) = {
      if (false) Util.time("Exact diagonalization") {
        val eig = KPM.eigenvaluesExact(q.matrix)
        val action  = eig.filter(_ < mup).map(_ - mup).sum
        val filling = eig.filter(_ < mup).size.toDouble / eig.size
        (eig, Array[Float](), action, filling)
      }
      else Util.time("Exact moments") {
        val moments = kpm.momentsExact(order_exact)
        val action  = moments dot c2_action
        val filling = moments dot c2_filling
        (Array[Double](), moments, action, filling)
      }
    }
    println("  Action  = %.7g".format(action))
    println("  Filling = %g".format(filling))
    
    // dump configuration
    val time = iter*dumpPeriod*dt
    val snap = KondoSnap(time=time, action=action, filling=filling, eig=eig, spin=q.field, moments=moments)
    val fn = dumpFilename(iter)
    println("Dumping %s (t=%g)".format(fn, time))
    kip.util.Util.writeStringToFile(serialize(snap), fn)
  }
  
  // the energy scales used in the langevin equation are reduced by q.e_scale (typically ~10)
  // this effectively reduces "time" scale by q.e_scale
  // therefore, to get the actual parameter (z=dt/nrand), one must divide by q.e_scale 
  val lang = new OverdampedLangevin(x=q.field, T=Tp, dt=dt, subIter=200, rand=kpm.rand) {
    override def calcForce(x: Array[R], f: Array[R]) {
      q.fillMatrix(q.matrix)
      val r = kpm.randomVector()
      kpm.functionAndGradient(r, c_action, q.delMatrix)
      q.fieldDerivative(q.delMatrix, f)
    }
    override def projectToBase(x: Array[R]) {
      q.normalizeField(x, validate=true)
    }
  }
  
  if (iter == 0) {
    calcMomentsAndDump()
    iter += 1
  }
  
  while (true) {
    Util.time("Langevin dynamics") (for (_ <- 0 until dumpPeriod) {
      lang.step()
      val hermitDev = math.sqrt((q.matrix - q.matrix.dag).norm2.abs)
      require(hermitDev < 1e-6, "Found non-hermitian hamiltonian! Deviation: "+hermitDev)
    })
    calcMomentsAndDump()
    iter += 1
  }
}
