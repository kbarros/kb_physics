package kip.projects.quantum

import net.liftweb.json
import java.io.File
import kip.enrich._
import kip.math.Vec3


object KondoStats extends App {
  val dir = args(0)
  val conf = KondoApp.readConfig(new File(dir+"/cfg.json"))
  val w = conf.w
  val h = conf.h
  val dumpdir = new java.io.File(dir+"/dump")
  require(dumpdir.isDirectory(), "Cannot load directory %s".format(dumpdir))

  val time   = new collection.mutable.ArrayBuffer[Double]()
  val energy = new collection.mutable.ArrayBuffer[Double]()
  val chiral = new collection.mutable.ArrayBuffer[Double]()
  val a1 = new collection.mutable.ArrayBuffer[Double]()
  val a2 = new collection.mutable.ArrayBuffer[Double]()
  val a3 = new collection.mutable.ArrayBuffer[Double]()

  def readSpin(x: Int, y: Int, field: Array[R]): Vec3 = {
    require(x < w && y < h)
    val sx = field(0 + 3*(x + w*y))
    val sy = field(1 + 3*(x + w*y))
    val sz = field(2 + 3*(x + w*y))
    Vec3(sx, sy, sz)
  }

  def calcChirality(field: Array[R]): Double = {
    val stats = new kip.util.Statistics.OnlineVariance
    for (y <- (h-2) to 0 by -1;
         x <- 0 until w-1) {
      //
      // s3 - s4
      //  \ /  \
      //  s1 - s2
      //
      val s1 = readSpin((x+0)%w, (y+0)%h, field)
      val s2 = readSpin((x+1)%w, (y+0)%h, field)
      val s3 = readSpin((x+0)%w, (y+1)%h, field)
      val s4 = readSpin((x+1)%w, (y+1)%h, field)

      stats.accum(s1 dot (s4 cross s3))
      stats.accum(s1 dot (s2 cross s4))
    }
    stats.mean
  }
  
  def calcCorrelations(field: Array[R]): (Double, Double, Double) = {
    var a1 = 0.0
    var a2 = 0.0
    var a3 = 0.0
    for (y <- 0 until h;
         x <- 0 until w) {
      //
      // s3 - s4
      //  \ /  \
      //  s1 - s2
      //
      val s1 = readSpin((x+0)%w, (y+0)%h, field)
      val s2 = readSpin((x+1)%w, (y+0)%h, field)
      val s3 = readSpin((x+0)%w, (y+1)%h, field)
      val s4 = readSpin((x+1)%w, (y+1)%h, field)
      a1 += s1 dot s2 
      a2 += s1 dot s4
      a3 += s1 dot s3
    }
    (a1/(w*h), a2/(w*h), a3/(w*h))
  }
  
  for (f <- dumpdir.listFiles() /* ; if i < 50 */ ) {
    implicit val formats = json.DefaultFormats
    val snap = json.Serialization.read[KondoSnap](f.slurp)
    time += snap.time
    energy += snap.action
    chiral += calcChirality(snap.spin)
    
    val as = calcCorrelations(snap.spin)
    a1 += as._1
    a2 += as._2
    a3 += as._3
  }
  
  val i = 0 // time.indexWhere(_ > (time.last / 5)) // skip first 1/5=20%
  val data = new scikit.dataset.PointSet(time.drop(i).toArray, chiral.drop(i).toArray)
  scikit.util.Commands.plot(data)
  
  val ba1 = new kip.util.BlockAnalysis(energy.drop(i).toArray)
  val ba2 = new kip.util.BlockAnalysis(chiral.drop(i).toArray)
  
  println("T energy +- chirality +- %n%f %.8g %.8g %.8g %.8g".format(conf.T, ba1.mean, ba1.error, ba2.mean, ba2.error)) // conf.dt_per_rand/10
  
  scikit.util.Commands.plot(new scikit.dataset.PointSet(time.toArray, a1.toArray))
  scikit.util.Commands.replot(new scikit.dataset.PointSet(time.toArray, a2.toArray))
  scikit.util.Commands.replot(new scikit.dataset.PointSet(time.toArray, a3.toArray))

  println("%g %g %g".format(a1.sum, a2.sum, a3.sum))
  
  if (!ba1.isDecorrelated) {
    println("not decorrelated!")
    ba1.blocks.foreach(b => println(b))
  }

  // println("normalize chirality by %g = 4 / 3^(3/2)".format(4 / math.pow(3, 1.5)))
}
