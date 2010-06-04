package kip.projects.dlc

import kip.util.{LammpsParser, Snapshot}
import kip.util.Util._

object Nano {
  // atom types 
  val typPatch = 1
  val typCore = 2
  val typSphereCation = 3
  val typCation = 4
  val typAnion = 5
  
  
  def averageTemperature(snaps: Seq[Snapshot]) = {
    snaps.map(_.thermo.temperature).sum / snaps.size
  }
  
  def pairCorrelationBins(dr: Double, rmax: Double) = {
    val nbins = (rmax/dr).toInt
    Array.range(0, nbins).map(bin => (bin + 0.5) * dr)
  }

  def pairCorrelation(snaps: Seq[Snapshot], dr: Double, rmax: Double, ids1: Seq[Int], ids2: Seq[Int]) = {
    val volume = snaps(0).volume
    val r = pairCorrelationBins(dr, rmax)
    val g = new Array[Double](r.size)
    
    // sum over all snapshots, and all pairs of particles. distances are binned into g
    for ((s,iter) <- snaps.zipWithIndex) {
      // if (iter % 100 == 0)
      //  println("Processing snapshot "+iter)
      for (i1 <- ids1; i2 <- ids2; if i1 != i2) {
        val dist = math.sqrt(s.distance2(i1, i2))
        val bin = (dist/dr).toInt
        if (bin < r.size)
          g(bin) += 1
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
    for (bin <- 0 until r.size) {
      val volume_fraction = 4*math.Pi*r(bin)*r(bin)*dr / volume
      g(bin) /= volume_fraction * snaps.size * uniquePairs
    }
    
    g
  }
  
  def pairCorrelationWithError(snaps: Seq[Snapshot], dr: Double, rmax: Double, ids1: Seq[Int], ids2: Seq[Int]) = {
    // partition snapshots into 199 equal sized groups
    val groupSize = math.max(snaps.size / 200, 1)
    val snapsGrouped = snaps.grouped(groupSize).toArray.dropRight(1)
    
    // for each distance 'r', construct array of approximations to g(r)
    val gs = snapsGrouped.map(pairCorrelation(_, dr, rmax, ids1, ids2)).transpose
    
    // estimate pair correlation as g_mean +- g_err
    gs.map (new kip.util.BlockAnalysis(_))
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
    
    val r = pairCorrelationBins(dr, rmax)
    
    // sphere-sphere correlation
    val b1 = {
      val s = snaps1(0)
      val idsCore = 0 until s.natoms
      time(pairCorrelationWithError(snaps1, dr, rmax, idsCore, idsCore), "Sphere-sphere")
    }
    
    // sphere-ion correlation
    val (b2, b3) = {
      val s = snaps2(0)
      def filterIds(f: Int => Boolean) = (0 until s.natoms) filter f
      val idsCore    = filterIds (i => s.typ(i) == typCore)
      val idsCation  = filterIds (i => s.typ(i) == typCation || s.typ(i) == typSphereCation)
      val idsAnion   = filterIds (i => s.typ(i) == typAnion)
      (time(pairCorrelationWithError(snaps2, dr, rmax, idsCore, idsCation), "Sphere-cation"),
       time(pairCorrelationWithError(snaps2, dr, rmax, idsCore, idsAnion), "Sphere-anion"))
    }
    
    if (b1.exists(b => b.error > 0 && !b.isDecorrelated))
      println("Sphere-sphere g(r) not decorrelated!")
    if (b2.exists(b => b.error > 0 && !b.isDecorrelated))
      println("Sphere-ion g(r) not decorrelated!")
    
    val formatted = formatDataInColumns(
      ("radii", r),
      ("g(core-core)", b1.map(_.mean)),
      ("err", b1.map(_.error)),
      ("g(core-cation)", b2.map(_.mean)),
      ("err", b2.map(_.error)),
      ("g(core-anion)", b3.map(_.mean)),
      ("err", b3.map(_.error)),
      ("cc_err_err", b1.map(_.error_error)),
      ("ci_err_err", b2.map(_.error_error))
    )
    writeStringToFile(formatted, "results.dat")
  }
}
