package kip.projects.dlc

import kip.util.{BlockAnalysis, LammpsParser, RangeArray, Snapshot}
import kip.util.Util._
import kip.util.Statistics.{mean}


object Two {
  // atom types 
  val typPatch = 1
  val typCore = 2
  val typSphereCation = 3
  val typCation = 4
  val typAnion = 5
  
  def averageTemperature(snaps: Seq[Snapshot]) = {
    snaps.foreach(s => if (s.thermo.temperature.isNaN) println("NaN temp at "+s.time))
    mean(snaps.map(_.thermo.temperature))
  }

  def averageEnergy(snaps: Seq[Snapshot]) = {
    snaps.foreach(s => if (s.thermo.energy.isNaN) println("NaN energy at "+s.time))
    mean(snaps.map(_.thermo.energy))
  }
  
  def accumulatedCharge(snaps: Seq[Snapshot], dx: Double, xmax: Double, ids: Seq[Int]) = {
    val volume = snaps(0).volume
    val ra = RangeArray.fill(-xmax, xmax, dx)(0.0)
    
    // sum over all snapshots, and all pairs of particles. distances are binned into g
    for ((s,iter) <- snaps.zipWithIndex) {
      val lx = s.hi.x - s.lo.x
      def wrap(x: Double) = if (x < lx/2) x else x - lx
      
      // if (iter % 100 == 0)
      //  println("Processing snapshot "+iter)
      for (i <- ids) {
        val x = wrap(s.x(i))
        if (ra.isDefinedAt(x)) {
          ra(x) += s.q(i) / snaps.size
        }
      }
    }
    
    ra
  }
  
  /*
  def accumulatedChargeWithError(snaps: Seq[Snapshot], dx: Double, xmax: Double, ids: Seq[Int]) = {
    // partition snapshots into 199 equal sized groups
    val groupSize = math.max(snaps.size / 200, 1)
    val snapsGrouped = snaps.grouped(groupSize).toArray.dropRight(1)
    
    // for each distance 'r', construct array of approximations to g(r)
    val gs = RangeArray.transpose(snapsGrouped.map(accumulatedCharge(_, dx, xmax, ids)))
    
    // estimate pair correlation as g_mean +- g_err
    gs.map (new kip.util.BlockAnalysis(_))
  }
  */
  
  def go(tbegin: Long, tmax: Long = Long.MaxValue, dx: Double, xmax: Double) {
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
    println("Average energy      = "+averageEnergy(snaps1))
    println("Processing "+snaps1.size+" snapshots")
    
    def filterIds(s: Snapshot, f: Int => Boolean) = (0 until s.natoms) filter f

    val cationCharge = {
      val s = snaps1(0)
      val idsCation  = filterIds (s, i => s.typ(i) == typSphereCation)
      time(accumulatedCharge(snaps1, dx, xmax, idsCation), "Cation")
    }

    val patchCharge = {
      val s = snaps2(0)
      val idsPatch = filterIds (s, i => s.typ(i) == typPatch)
      time(accumulatedCharge(snaps2, dx, xmax, idsPatch),  "Patch")
    }
    
    // sphere-ion correlation
    val formatted = formatDataInColumns(
      ("pos", patchCharge.elemCenters),
      ("cationCharge", cationCharge.elems),
      ("patchCharge",  patchCharge.elems)
    )
    writeStringToFile(formatted, "results.dat")
    println()
  }
}
