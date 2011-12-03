package kip.projects.quantum

import net.liftweb.json
import java.io.File
import kip.util.Util
import kip.enrich._
import ctor._

case class KondoConf(w: Int, h: Int, t: Double, J_H: Double, mu: Double,
                     order: Int, de: Double, dt_per_rand: Double, nrand: Int, dumpPeriod: Int,
                     initConf: String)
case class KondoSnap(time: Double, action: Double, filling: Double, eig: Array[Double], moments: Array[Float], spin: Array[Float])


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
  
  val q = new Quantum(w=w, h=h, t=t, J_H=J_H, e_min= -10, e_max= 10)

  val kpm = try {
    import kip.projects.cuda._
    new CuKPM(new JCudaWorld(deviceIndex), q.matrix, nrand)
  } catch {
    case _ => new KPM(q.matrix, nrand)
  }
  
  initConf match {
    case "random" => q.setFieldRandom(q.field, kpm.rand)
    case "allout" => q.setFieldAllOut(q.field)
    case "threeout" => q.setFieldThreeOut(q.field)
  }
  q.fillMatrix(q.matrix)
  val dt = dt_per_rand * nrand
  val mup = q.scaleEnergy(mu)
  
  val fn_action:  (R => R) = e => if (e < mup) (e - mup) else 0
  val fn_filling: (R => R) = e => if (e < mup) (1.0 / q.matrix.numRows) else 0
  val c = KPM.expansionCoefficients(order, de, fn_action)
  
  println("N=%d matrix, %d moments".format(q.matrix.numRows, order))
  
  for (iter <- 0 until 1000) {
    Util.time("Iteration "+iter) (for (iter2 <- 0 until dumpPeriod) {
      val r = kpm.randomVector()
      val f0 = kpm.functionAndGradient(r, c, q.delMatrix)
      q.fieldDerivative(q.delMatrix, q.delField)
      for (i <- q.field.indices) {
        q.field(i) -= dt * q.delField(i)
      }
      q.normalizeField(q.field, validate=true)
      q.fillMatrix(q.matrix)
      require(math.sqrt((q.matrix - q.matrix.dag).norm2.abs) < 1e-14, "Found non-hermitian hamiltonian!")
    })

    // exact moments
    val (eig, moments, action, filling) = {
      if (false) Util.time("Exact diagonalization") {
        val eig = KPM.eigenvaluesExact(q.matrix)
        val action  = eig.filter(_ < mup).map(_ - mup).sum
        val filling = eig.filter(_ < mup).size.toDouble / eig.size
        (eig, Array[Float](), action, filling)
      }
      else Util.time("Exact moments") {
        val order2 = 1024
        val moments = kpm.momentsExact(order2)
        val action  = moments dot KPM.expansionCoefficients(order2, de, fn_action)
        val filling = moments dot KPM.expansionCoefficients(order2, de, fn_filling)
        (Array[Double](), moments, action, filling)
      }
    }
    println("  Action  = %g".format(action))
    println("  Filling = %g".format(filling))

    val time = iter*dumpPeriod*dt
//    val snap = KondoSnap(time=time, action=action, filling=filling, eig=eig, spin=q.field)
    val snap = KondoSnap(time=time, action=action, filling=filling, eig=eig, spin=q.field, moments=moments)
    implicit val formats = json.DefaultFormats
    val serialized = json.Serialization.write(snap)
    kip.util.Util.writeStringToFile(serialized, dumpdir+"/%04d.json".format(iter))
  }
}
