package kip.lmps


object Nano {
  def main(args: Array[String]) {
    go("/Users/kbarros/dev/repo/projects/dielectric/nano.L15/dump.dat", 1000)
  }
  
  // types 
  val typPatch = 1
  val typCore = 2
  val typCation = 3
  val typAnion = 4
  

  def pairCorrelation(snaps: Seq[Snapshot], dr: Double, rmax: Double, typ1: Int, typ2: Int) = {
    val types = snaps(0).typ
    val volume = snaps(0).volume
    val indices1 = types.indices.filter(types(_) == typ1)
    val indices2 = types.indices.filter(types(_) == typ2)
    val nbins = (rmax/dr).toInt
    val g = new Array[Double](nbins)
    
    // sum over all snapshots, and all pairs of particles. distances are binned into g
    for ((s,iter) <- snaps.zipWithIndex) {
      if (iter % 100 == 0)
        println("Processing snapshot "+iter)
      for (i1 <- indices1; i2 <- indices2; if i1 != i2) {
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
      g(bin) /= volume_fraction * snaps.size * indices1.size * indices2.size
    }
    
    (r, g)
  }
  
  def go(fname: String, tbegin: Long) {
    val snaps = LammpsParser.readLammpsDump(fname) filter {_.time > tbegin}
    // LammpsParser.weaveThermoData(snaps, LammpsParser.readLammpsThermo("log.lammps"))
    
    println("Processing "+snaps.size+" snapshots")
    
    val (r, g1) = pairCorrelation(snaps, dr=0.2, rmax=12, typCore, typCore)
    val (_, g2) = pairCorrelation(snaps, dr=0.2, rmax=12, typCore, typCation)
    val (_, g3) = pairCorrelation(snaps, dr=0.2, rmax=12, typCation, typCation)
    
    val formatted = Util.formatDataInColumns(
      ("radii", r),
      ("g(core-core)", g1),
      ("g(core-cation)", g2),
      ("g(cation-cation)", g3)
    )
    Util.writeStringToFile(formatted, "/Users/kbarros/Desktop/results.dat")
    // print(formatted)
  }
}
