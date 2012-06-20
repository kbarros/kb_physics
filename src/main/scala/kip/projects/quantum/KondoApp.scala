package kip.projects.quantum

import net.liftweb.json
import java.io.File
import kip.util.Util
import kip.enrich._
import ctor._
import scala.math._

case class KondoConf(w: Int, h: Int, t: Double, J_H: Double, B_n: Int, T: Double, mu: Double,
                     order: Int, order_exact: Int, de: Double, dt_per_rand: Double,
                     nrand: Int, dumpPeriod: Int, initConf: String)
case class KondoSnap(time: Double, action: Double, filling: Double, eig: Array[Double],
                     moments: Array[Float], spin: Array[Float])


object KondoApp extends App {
  def readConfig(f: File): KondoConf = {
    val confStripped = f.slurp.split("\n") map { _.split("""//""")(0) } mkString("\n")
    implicit val formats = json.DefaultFormats
    json.Serialization.read[KondoConf](confStripped)
  }
  
  def dumpFile(iter: Int) = dumpdir+"/%04d.json".format(iter)
  
  if (args.size != 2) {
    println("KondoApp requires <dir> and <device> parameters")
    sys.exit
  }
  val dir = args(0)
  val deviceIndex = args(1).toInt
  
  val conf = readConfig(new File(dir+"/cfg.json"))
  import conf._
  
  // create output directory for spin configurations
  val dumpdir = new java.io.File(dir+"/dump")
  Util.createEmptyDir(dumpdir)
  
  val q = new Quantum(w=w, h=h, t=t, J_H=J_H, B_n=B_n, e_min= -10, e_max= 10)

  val kpm = try {
    import kip.projects.cuda._
    new CuKPM(new JCudaWorld(deviceIndex), q.matrix, nrand)
  } catch {
    case _ => {
      println("CUDA device not available! Defaulting to CPU implementation.")
      new KPM(q.matrix, nrand)
    }
  }
  
  val iter0 = initConf match {
    case "random" => {q.setFieldRandom(q.field, kpm.rand); 0}
    case "ferro" => {q.setFieldFerro(q.field); 0}
    case "allout" => {q.setFieldAllOut(q.field); 0}
    case "threeout" => {q.setFieldThreeOut(q.field); 0}
    case s => {
      val lastDump = s.toInt
      val f = new File(dumpFile(lastDump))
      println("Reading initial configuration from file: "+f)
      implicit val formats = json.DefaultFormats
      val snap = json.Serialization.read[KondoSnap](f.slurp)
      snap.spin.copyToArray(q.field)
      lastDump+1
    }
  }
  
  // don't overwrite any existing data
  require(!(new File(dumpFile(iter0))).exists, "Refuse to overwrite dump file %s".format(dumpFile(iter0)))
  
  q.fillMatrix(q.matrix)
  val dt = dt_per_rand * nrand
  val mup = q.scaleEnergy(mu)
  val Tp = q.scaleEnergy(T)
  
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
  
  val c_action   = KPM.expansionCoefficients(order, de, fn_action)
  val c2_action  = KPM.expansionCoefficients(order_exact, de, fn_action)
  val c2_filling = KPM.expansionCoefficients(order_exact, de, fn_filling)
  
  println("N=%d matrix, %d moments".format(q.matrix.numRows, order))
  
  val lang = new OverdampedLangevin(x=q.field, T=Tp, dt=dt, subIter=10, rand=kpm.rand) {
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
  
  for (iter <- iter0 until 1000) {
    Util.time("Iteration "+iter) (for (iter2 <- 0 until dumpPeriod) {
      lang.step()
      val hermitDev = math.sqrt((q.matrix - q.matrix.dag).norm2.abs)
      require(hermitDev < 1e-6, "Found non-hermitian hamiltonian! Deviation: "+hermitDev)
    })

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

    val time = (iter+1)*dumpPeriod*dt
//    val snap = KondoSnap(time=time, action=action, filling=filling, eig=eig, spin=q.field)
    val snap = KondoSnap(time=time, action=action, filling=filling, eig=eig, spin=q.field, moments=moments)
    implicit val formats = json.DefaultFormats
    val serialized = json.Serialization.write(snap)
    kip.util.Util.writeStringToFile(serialized, dumpFile(iter))
  }
}
