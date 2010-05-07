package kip.lmps


object Nano {
  def main(args: Array[String]) {
    go("/Users/kbarros/dev/repo/projects/dielectric/nano.L15/dump.dat", 2000, 0.1, 10)
  }
  
  // atom types 
  val typPatch = 1
  val typCore = 2
  val typCation = 3
  val typAnion = 4
  
  
  def averageTemperature(snaps: Seq[Snapshot]) = {
    snaps.map(_.thermo.temperature).sum / snaps.size
  }

  def pairCorrelation(snaps: Seq[Snapshot], dr: Double, rmax: Double, ids1: Seq[Int], ids2: Seq[Int]) = {
    val volume = snaps(0).volume
    val nbins = (rmax/dr).toInt
    val g = new Array[Double](nbins)
    
    // sum over all snapshots, and all pairs of particles. distances are binned into g
    for ((s,iter) <- snaps.zipWithIndex) {
      if (iter % 100 == 0)
        println("Processing snapshot "+iter)
      for (i1 <- ids1; i2 <- ids2; if i1 != i2) {
        val dist = math.sqrt(s.distance2(i1, i2))
        val bin = (dist/dr).toInt
        if (bin < nbins)
          g(bin) += 1
      }
    }

    // bin radii (average)
    val r = Array.range(0, nbins).map(bin => (bin + 0.5) * dr)
    
    // normalize pair-correlation to unity for homogeneous density
    for (bin <- 0 until nbins) {
      val volume_fraction = 4*math.Pi*r(bin)*r(bin)*dr / volume
      g(bin) /= volume_fraction * snaps.size * ids1.size * ids2.size
    }
    
    (r, g)
  }
  
  def go(fname: String, tbegin: Long, dr: Double, rmax: Double) {
    val snaps = LammpsParser.readLammpsDump(fname) filter {_.time > tbegin}
    LammpsParser.weaveThermoData(snaps, LammpsParser.readLammpsThermo("log.lammps"))
    println("Average temperature = "+averageTemperature(snaps))
    println("Processing "+snaps.size+" snapshots")
    
    val s = snaps(0)
    val types = snaps(0).typ
    def isPositive(q: Array[Double], i: Int) = if (q == null) true else q(i) > 0
    def isNegative(q: Array[Double], i: Int) = if (q == null) true else q(i) < 0
    def filterIds(f: Int => Boolean) = (0 until s.natoms) filter f
    
    val idsCorePos = filterIds (i => s.typ(i) == typCore && isPositive(s.q, i))
    val idsCoreNeg = filterIds (i => s.typ(i) == typCore && isNegative(s.q, i))
    val idsCore    = filterIds (i => s.typ(i) == typCore)
    val idsCation  = filterIds (i => s.typ(i) == typCation)
    val idsAnion   = filterIds (i => s.typ(i) == typAnion)
    
    val (r, g1) = pairCorrelation(snaps, dr, rmax, idsCorePos, idsCoreNeg)
    val (_, g2) = pairCorrelation(snaps, dr, rmax, idsCore, idsCation)
    val (_, g3) = pairCorrelation(snaps, dr, rmax, idsCation, idsCation)
    
    val formatted = Util.formatDataInColumns(
      ("radii", r),
      ("g(core-core)", g1),
      ("g(core-cation)", g2),
      ("g(cation-cation)", g3)
    )
    Util.writeStringToFile(formatted, "results.dat")
    // print(formatted)
  }
}
