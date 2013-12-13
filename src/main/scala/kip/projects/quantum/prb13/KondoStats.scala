package kip.projects.quantum.prb13

import java.io.File
import scala.Float.float2double

import kip.projects.quantum._
import kip.enrich._
import kip.math.Vec3
import kip.util.Statistics._

// For generation of phase diagram in PRB

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
    val snap = kip.util.JacksonWrapper.deserialize[KondoSnap](f.slurp)
    time += snap.time
    energy += snap.action
    chiral += calcChirality(snap.spin)
    
    val as = calcCorrelations(snap.spin)
    a1 += as._1
    a2 += as._2
    a3 += as._3
  }
  
  val i = time.indexWhere(_ > (time.last / 4)) // skip first 1/4=25%
  val data = new scikit.dataset.PointSet(time.drop(i).toArray, chiral.drop(i).toArray)
  scikit.util.Commands.plot(data)
  
  val ba1 = new kip.util.BlockAnalysis(energy.drop(i).toArray)
  val ba2 = new kip.util.BlockAnalysis(chiral.drop(i).toArray)
  
  // printf("T energy +- chirality +-\n%f %.8g %.8g %.8g %.8g\n\n",conf.T, ba1.mean, ba1.error, ba2.mean, ba2.error)
  printf("T energy +- chirality +-\n%f %.8g %.8g %.8g %.8g\n\n",conf.T, ba1.mean, ba1.error, mean(chiral.drop(i)),stddev(chiral.drop(i)))
  
  val groupSize = 1 // 13 // 50
  val tFilter = 2e5 // 4e6 
  
  val a1_avg = a1.grouped(groupSize).map(mean(_)).toArray
  val a2_avg = a2.grouped(groupSize).map(mean(_)).toArray
  val a3_avg = a3.grouped(groupSize).map(mean(_)).toArray
  val t_avg = time.grouped(groupSize).map(mean(_)).toArray
  
  val a1_stddev = a1.grouped(groupSize).map(stddev(_)).toArray
  val a2_stddev = a2.grouped(groupSize).map(stddev(_)).toArray
  val a3_stddev = a3.grouped(groupSize).map(stddev(_)).toArray
  
  
  val (b1, b2, b3) = (for (i <- t_avg.indices;
                      if t_avg(i) > tFilter) yield {
    val vs = List((a1_avg(i), a1_stddev(i)), (a2_avg(i), a2_stddev(i)), (a3_avg(i), a3_stddev(i)))
    val vs_sorted = vs.sortWith { case ((avg1, _), (avg2, _)) => avg1 < avg2 }
    val List(b1, b2, b3) = vs_sorted
    (b1, b2, b3)
  }).unzip3
  
  val t_avg2 = t_avg.filter(_ > tFilter) 
    
  scikit.util.Commands.plot(new scikit.dataset.PointSet(time.toArray, a1.toArray))
  scikit.util.Commands.replot(new scikit.dataset.PointSet(time.toArray, a2.toArray))
  scikit.util.Commands.replot(new scikit.dataset.PointSet(time.toArray, a3.toArray))
  
  scikit.util.Commands.plot(new scikit.dataset.PointSet(t_avg2, b1.map(_._1).toArray))
  scikit.util.Commands.replot(new scikit.dataset.PointSet(t_avg2, b2.map(_._1).toArray))
  scikit.util.Commands.replot(new scikit.dataset.PointSet(t_avg2, b3.map(_._1).toArray))
  
  printf("T b1 +- b2 +- b3 +-\n")
  printf("%f %g %g %g %g %g %g",
      conf.T,
      mean(b1.map(_._1)), mean(b1.map(_._2)),
      mean(b2.map(_._1)), mean(b2.map(_._2)),
      mean(b3.map(_._1)), mean(b3.map(_._2)))
//  printf(" # tFilter=%g, groupSize=%dx%g\n", tFilter, groupSize, time(1)-time(0))
  
//  if (!ba1.isDecorrelated) {
//    println("not decorrelated!")
//    ba1.blocks.foreach(b => println(b))
//  }

  // println("normalize chirality by %g = 4 / 3^(3/2)".format(4 / math.pow(3, 1.5)))
}
