package kip.projects.quantum

import net.liftweb.json
import java.io.File
import kip.enrich._


object KondoStats extends App {
  val dir = args(0)
  val conf = KondoApp.readConfig(new File(dir+"/cfg.json"))
  val w = conf.w
  val h = conf.h
  val dumpdir = new java.io.File(dir+"/dump")
  require(dumpdir.isDirectory(), "Cannot load directory %s".format(dumpdir))

  val time   = new collection.mutable.ArrayBuffer[Double]()
  val energy = new collection.mutable.ArrayBuffer[Double]()
  
  for (f <- dumpdir.listFiles() /* ; if i < 50 */ ) {
    implicit val formats = json.DefaultFormats
    val snap = json.Serialization.read[KondoSnap](f.slurp)
    time += snap.time
    energy += snap.action
  }
  
  
  val i = 0 // time.indexWhere(_ > 5000)
  val data = new scikit.dataset.PointSet(time.drop(i).toArray, energy.drop(i).toArray)
  scikit.util.Commands.plot(data)
  
  val ba = new kip.util.BlockAnalysis(energy.drop(i).toArray)
  println("energy = %.8g +- %.8g".format(ba.mean, ba.error))
  if (!ba.isDecorrelated) {
    println("not decorrelated!")
    ba.blocks.foreach(b => println(b))
  }
}
