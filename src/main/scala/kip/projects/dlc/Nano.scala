package kip.projects.dlc

import kip.util.{LammpsParser, Snapshot, BlockAnalysis}
import kip.util.{RangeArray}
import kip.util.Util._
import kip.math.Vec3


object Nano {
  // atom types 
  val typPatch = Seq(1)
  val typCore = Seq(2)
  val typCation = Seq(3, 4) // include "sphere" and "salt" cations
  val typAnion = Seq(5)
  
  def averageTemperature(snaps: Seq[Snapshot]) = {
    snaps.foreach(s => if (s.thermo.temperature.isNaN) println("NaN temp at "+s.time))
    snaps.map(_.thermo.temperature).sum / snaps.size
  }
  
  def bondAngle(dx1: Vec3, dx2: Vec3): Double = {
    math.acos((dx1 dot dx2) / (dx1.norm*dx2.norm))
  }

  def sizeTwoSubsets[A](a: Seq[A]): Seq[(A,A)] = {
    for (i <- 0   until a.size;
         j <- i+1 until a.size) yield {
      (a(i), a(j))
    }
  }
  
  def bondAngleHistogram(snaps: Seq[Snapshot], dtheta: Double, typs1: Seq[Int], typs2: Seq[Int], rcutoff: Double) = {
    val g = RangeArray.fill(xmin=0, xmax=math.Pi, dx=dtheta)(0d)
    for ((s,iter) <- snaps.zipWithIndex) {
      val ids1 = (0 until s.natoms) filter (i => typs1.contains(s.typ(i)))
      val ids2 = (0 until s.natoms) filter (i => typs2.contains(s.typ(i)))
      for (i <- ids1) {
        val bonded = ids2.filter { j => (i != j) && (s.distance2(i, j) < rcutoff*rcutoff) }
        for ((j1, j2) <- sizeTwoSubsets(bonded)) {
          val d1 = s.displacement(j1, i)
          val d2 = s.displacement(j2, i)
          g(bondAngle(d1, d2)) += 1.0
        }
      }
    }

    // normalize
    val ids1 = (0 until snaps(0).natoms) filter (i => typs1.contains(snaps(0).typ(i)))
    for (j <- g.elemCenters) {
      g(j) /= dtheta * snaps.size * ids1.size * (2*math.Pi*math.sin(j))
    }
    
    g
  }
  
  def pairCorrelation(snaps: Seq[Snapshot], dr: Double, rmax: Double, typs1: Seq[Int], typs2: Seq[Int]): RangeArray[Double] = {
    val volume = snaps(0).volume
    val g = RangeArray.fill(xmin=0, xmax=rmax, dx=dr)(0d)

    // sum over all snapshots, and all pairs of particles. distances are binned into g
    for ((s,iter) <- snaps.zipWithIndex) {
      val ids1 = (0 until s.natoms) filter (i => typs1.contains(s.typ(i)))
      val ids2 = (0 until s.natoms) filter (i => typs2.contains(s.typ(i)))
      // if (iter % 100 == 0)
      //  println("Processing snapshot "+iter)
      for (i1 <- ids1; i2 <- ids2; if i1 != i2) {
        val r = math.sqrt(s.distance2(i1, i2))
        if (g.isDefinedAt(r))
          g(r) += 1
      }
      
      // // while loops give no apparent speedup
      // var i1 = 0
      // while (i1 < ids1.size) {
      //   var i2 = 0
      //   while (i2 < ids2.size) {
      //     if (i1 != i2) {
      //       val dist = math.sqrt(s.distance2(ids1(i1), ids2(i2)))
      //       val bin = (dist/dr).toInt
      //       if (bin < r.size)
      //         g(bin) += 1
      //     }
      //     i2 += 1
      //   }
      //   i1 += 1
      // }
    }

    // normalize pair-correlation so that it becomes unity at homogeneous density
    val ids1 = (0 until snaps(0).natoms) filter (i => typs1.contains(snaps(0).typ(i)))
    val ids2 = (0 until snaps(0).natoms) filter (i => typs2.contains(snaps(0).typ(i)))
    val uniquePairs = (ids1.size * ids2.size) - (ids1.toSet & ids2.toSet).size
    for (r <- g.elemCenters) {
      val volume_fraction = 4*math.Pi*r*r*dr / volume
      g(r) /= volume_fraction * snaps.size * uniquePairs
    }
    
    g
  }
  
  def pairCorrelationWithError(snaps: Seq[Snapshot], dr: Double, rmax: Double, typs1: Seq[Int], typs2: Seq[Int]): RangeArray[BlockAnalysis] = {
    // partition snapshots into 199 equal sized groups
    val groupSize = math.max(snaps.size / 200, 1)
    
    // list of groups of snapshots
    val snapsGrouped = snaps.grouped(groupSize).toArray.dropRight(1)
    
    // list of approximations to function g(r)
    val gs = snapsGrouped.map(pairCorrelation(_, dr, rmax, typs1, typs2))

    // pair correlation estimate, g_mean +- g_err
    RangeArray.transpose(gs).map (new kip.util.BlockAnalysis(_))
  }
  
  
  def writeCorrelationFunctions(snaps2: Seq[Snapshot], dr: Double, rmax: Double) {
    val (b1, b2, b3, b4) = {
      val s0 = snaps2(0)
      (time("Sphere-sphere")(pairCorrelationWithError(snaps2, dr, rmax, typCore,   typCore)),
      time("Sphere-cation") (pairCorrelationWithError(snaps2, dr, rmax, typCore,   typCation)),
      time("Sphere-anion")  (pairCorrelationWithError(snaps2, dr, rmax, typCore,   typAnion)),
      time("Cation-cation") (pairCorrelationWithError(snaps2, dr, rmax, typCation, typCation)))
    }
    
    if (b1.exists(b => b.error > 0 && !b.isDecorrelated))
      println("Sphere-sphere g(r) not decorrelated!")
    if (b2.exists(b => b.error > 0 && !b.isDecorrelated))
      println("Sphere-ion g(r) not decorrelated!")
    
    def adjustedError(b: BlockAnalysis) =
      if (b.isDecorrelated) b.error else 5*b.error
    
    val formatted = formatDataInColumns(
      ("radii", b1.elemCenters),
      ("g(core-core)", b1.map(_.mean).elems),
      ("err", b1.map(adjustedError _).elems),
      ("g(core-cation)", b2.map(_.mean).elems),
      ("err", b2.map(adjustedError _).elems),
      ("g(core-anion)", b3.map(_.mean).elems),
      ("err", b3.map(adjustedError _).elems),
      ("g(cation-cation)", b4.map(_.mean).elems),
      ("err", b4.map(adjustedError _).elems)
    )
    writeStringToFile(formatted, "results.dat")
    println()
  }

  
  def writeAngleHistogram(snaps: Seq[Snapshot], dtheta: Double) {
    val s0 = snaps(0)
    val g = time("Angle histogram")(bondAngleHistogram(snaps, dtheta, typCore, typCore, rcutoff=8.5))
    val formatted = formatDataInColumns(
      ("radii", g.elemCenters),
      ("g(theta)", g.elems)
    )
    writeStringToFile(formatted, "angles.dat")
    println()
  }


  def go(tmin:Int, tmax:Int, dr: Double, rmax: Double, dtheta: Double, readEvery: Int) {
    // discard unused arrays to save space; filter by time
    def process(s: Snapshot) = {
      if (s.time > tmin) {
        s.id = null
        s.ix = null
        s.iy = null
        s.iz = null
        s.vx = null
        s.vy = null
        s.vz = null
        Some(s)
      }
      else {
        None
      }
    }
    def terminate(snaps: Seq[Snapshot]): Boolean = {
      snaps.lastOption.map(_.time > tmax).getOrElse(false)
    }
    // val snaps1 = time("Reading dump1.gz")(LammpsParser.readLammpsDump("dump1-0.gz", process, terminate, readEvery))
    val snaps2 = time("Reading dump2.gz")(LammpsParser.readLammpsDump("dump2-0.gz", process, terminate, readEvery))
    time("Weaving thermo")(LammpsParser.weaveThermoData(snaps2, LammpsParser.readLammpsThermo("log.lammps")))
    println("Processing "+snaps2.size+" of "+(snaps2.size*readEvery)+" snapshots")
    println("Average temperature = "+averageTemperature(snaps2))
    
    writeCorrelationFunctions(snaps2, dr, rmax)
    writeAngleHistogram(snaps2, dtheta=dtheta)
  }
  
  def main(args: Array[String]) {
    import kip.util.JacksonWrapper._
    
    val files = {
      if (args.size == 0)
        List("cfg.json", "../cfg.json", "../../cfg.json")
      else 
        args.toList
    }.map(new java.io.File(_))
    
    var params = Map[Any, Any]()
    for (file <- files) {
      if (file.exists) {
        println("Loading configuration file '"+file+"'")
        params ++= deserialize[Map[Any,Any]](kip.util.Util.readStringFromFile(file.toString))
      }
    }
    if (params.isEmpty) {
      println("Empty configuration")
      System.exit(1)
    }

    go(tmin=params("tmin").asInstanceOf[Int],
       tmax=params("tmax").asInstanceOf[Int],
       dr=params("dr").asInstanceOf[Double],
       rmax=params("rmax").asInstanceOf[Double],
       dtheta=params("dtheta").asInstanceOf[Double],
       readEvery=params("readEvery").asInstanceOf[Int])
  }
}
