package kip.projects.dlc

import kip.util.{LammpsParser, Snapshot, BlockAnalysis}
import kip.util.{RangeArray}
import kip.util.Util._
import kip.math.Vec3


object Nano {
  // atom types 
  val typPatch = 1
  val typCore = 2
  val typSphereCation = 3
  val typCation = 4
  val typAnion = 5
  
  def averageTemperature(snaps: Seq[Snapshot]) = {
    snaps.foreach(s => if (s.thermo.temperature.isNaN) println("NaN temp at "+s.time))
    snaps.map(_.thermo.temperature).sum / snaps.size
  }
  
  def bondAngle(dx1: Vec3, dx2: Vec3): Double = {
    math.acos((dx1 dot dx2) / (dx1.norm*dx2.norm))
  }

  def bondAngleHistogram(snaps: Seq[Snapshot], dj: Double, jmax: Double, ids1: Seq[Int], ids2: Seq[Int], rcutoff: Double) = {
    val g = RangeArray.fill(xmin=0, xmax=jmax, dx=dj)(0d)
    for ((s,iter) <- snaps.zipWithIndex) {
      for (i <- ids1) {
        val bonded = ids2.filter { j => (i != j) && (s.distance2(i, j) < rcutoff*rcutoff) }
        for (j1 <- 0 until bonded.size;
             j2 <- j1+1 until bonded.size) {
          val d1 = s.displacement(j1, i)
          val d2 = s.displacement(j2, i)
          g(bondAngle(d1, d2)) += 1.0
        }
      }
    }

    // normalize
    for (j <- g.binCenters) {
      g(j) /= dj * snaps.size * ids1.size
    }
    
    g
  }
  
  def pairCorrelation(snaps: Seq[Snapshot], dr: Double, rmax: Double, ids1: Seq[Int], ids2: Seq[Int]): RangeArray[Double] = {
    val volume = snaps(0).volume
    val g = RangeArray.fill(xmin=0, xmax=rmax, dx=dr)(0d)
    
    // sum over all snapshots, and all pairs of particles. distances are binned into g
    for ((s,iter) <- snaps.zipWithIndex) {
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
    val uniquePairs = (ids1.size * ids2.size) - (ids1.toSet & ids2.toSet).size
    for (r <- g.binCenters) {
      val volume_fraction = 4*math.Pi*r*r*dr / volume
      g(r) /= volume_fraction * snaps.size * uniquePairs
    }
    
    g
  }
  
  def pairCorrelationWithError(snaps: Seq[Snapshot], dr: Double, rmax: Double, ids1: Seq[Int], ids2: Seq[Int]): RangeArray[BlockAnalysis] = {
    // partition snapshots into 199 equal sized groups
    val groupSize = math.max(snaps.size / 200, 1)
    
    // list of groups of snapshots
    val snapsGrouped = snaps.grouped(groupSize).toArray.dropRight(1)
    
    // list of approximations to function g(r)
    val gs = snapsGrouped.map(pairCorrelation(_, dr, rmax, ids1, ids2))
    
    // pair correlation estimate, g_mean +- g_err
    RangeArray.transpose(gs).map (new kip.util.BlockAnalysis(_))
  }
  
  
  def go(tbegin: Long, tmax: Long = Long.MaxValue, dr: Double, rmax: Double) {
    // discard unused arrays to save space; filter by time
    def process(s: Snapshot) = {
      if (s.time > tbegin && s.time < tmax) {
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
    val snaps1 = time(LammpsParser.readLammpsDump("dump1.gz", process), "Reading dump1.gz")
    val snaps2 = time(LammpsParser.readLammpsDump("dump2.gz", process), "Reading dump2.gz")
    time(LammpsParser.weaveThermoData(snaps1, LammpsParser.readLammpsThermo("log.lammps")), "Weaving thermo")
    println("Average temperature = "+averageTemperature(snaps1))
    println("Processing "+snaps1.size+" snapshots")
    
    // sphere-sphere correlation
    val b1 = {
      val s = snaps1(0)
      val idsCore = 0 until s.natoms
      time(pairCorrelationWithError(snaps1, dr, rmax, idsCore, idsCore), "Sphere-sphere")
    }
    
    // sphere-ion correlation
    val (b2, b3, b4) = {
      val s = snaps2(0)
      def filterIds(f: Int => Boolean) = (0 until s.natoms) filter f
      val idsCore    = filterIds (i => s.typ(i) == typCore)
      val idsCation  = filterIds (i => s.typ(i) == typCation || s.typ(i) == typSphereCation)
      val idsAnion   = filterIds (i => s.typ(i) == typAnion)
      (time(pairCorrelationWithError(snaps2, dr, rmax, idsCore, idsCation), "Sphere-cation"),
       time(pairCorrelationWithError(snaps2, dr, rmax, idsCore, idsAnion), "Sphere-anion"),
       time(pairCorrelationWithError(snaps2, dr, rmax, idsCation, idsCation), "Cation-cation"))
    }
    
    if (b1.exists(b => b.error > 0 && !b.isDecorrelated))
      println("Sphere-sphere g(r) not decorrelated!")
    if (b2.exists(b => b.error > 0 && !b.isDecorrelated))
      println("Sphere-ion g(r) not decorrelated!")
    
    def adjustedError(b: BlockAnalysis) =
      if (b.isDecorrelated) b.error else 5*b.error
    
    val formatted = formatDataInColumns(
      ("radii", b1.binCenters),
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
}
