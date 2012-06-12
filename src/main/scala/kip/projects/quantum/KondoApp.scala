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
  
  initConf match {
    case "random" => q.setFieldRandom(q.field, kpm.rand)
    case "ferro" => q.setFieldFerro(q.field)
    case "allout" => q.setFieldAllOut(q.field)
    case "threeout" => q.setFieldThreeOut(q.field)
    case s => {
      val f = new File(s)
      println("Loading dump file: "+s)
      implicit val formats = json.DefaultFormats
      val snap = json.Serialization.read[KondoSnap](f.slurp)
    }
  }
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
  
  val lang = new OverdampedLangevin(x=q.field, T=Tp, dt=dt, subIter=1, rand=kpm.rand) {
    override def calcForce(x: Array[R], f: Array[R]) {
      q.fillMatrix(q.matrix)
      val r = kpm.randomVector()
      
      if (true)
        kpm.functionAndGradient(r, c_action, q.delMatrix) // approximate gradient
      else
        kpm.gradientExactDense(c_action, q.delMatrix) // exact gradient
      
      q.fieldDerivative(q.delMatrix, f)
    }
    override def projectToBase(x: Array[R], xp: Array[R]) {
      q.normalizeField(q.field, validate=true)
    }
  }
  
  for (iter <- 0 until 1000) {
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

    val time = iter*dumpPeriod*dt
//    val snap = KondoSnap(time=time, action=action, filling=filling, eig=eig, spin=q.field)
    val snap = KondoSnap(time=time, action=action, filling=filling, eig=eig, spin=q.field, moments=moments)
    implicit val formats = json.DefaultFormats
    val serialized = json.Serialization.write(snap)
    kip.util.Util.writeStringToFile(serialized, dumpdir+"/%04d.json".format(iter))
  }
}
