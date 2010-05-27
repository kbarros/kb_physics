package kip.projects.peb

import kip.util.{Snapshot, LammpsParser}
import kip.util.Vec3
import kip.util.Util.{sqr, formatDataInColumns, writeStringToFile}
import kip.util.Statistics.{mean}
import scala.math.sqrt


object PEB {
  def molecules(chainLength: Int, numChains: Int): Seq[Seq[Int]] = {
    for (i <- 0 until numChains)
      yield (i*chainLength) until ((i+1)*chainLength)
  }

  // unfortunately the image indices in much of the polymer brush data is invalid.
  // this function corrects the errors, making the polymers continuous
  def fixPolymerContinuity(snap: Snapshot, mols: Seq[Seq[Int]]) {
    val Lx = snap.hi.x - snap.lo.x
    val Ly = snap.hi.y - snap.lo.y
    val Lz = snap.hi.z - snap.lo.z
    
    for (poly <- mols) {
      for (i <- 0 until poly.size - 1) {
        val m1 = poly(i)
        val m2 = poly(i+1)
        
        def adjust(x1: Double, x2: Double, L: Double) = {
          var x2p = x2
          while (x2p - x1 > L/2)
            x2p -= L
          while (x2p - x1 < -L/2)
            x2p += L
          x2p
        }
        snap.x(m2) = adjust(snap.x(m1), snap.x(m2), Lx)
        snap.y(m2) = adjust(snap.y(m1), snap.y(m2), Ly)
        snap.z(m2) = adjust(snap.z(m1), snap.z(m2), Lz)
        
        // make sure each bond is small
        val d = sqrt(snap.getPoint(m1) distance2 snap.getPoint(m2))
        if (d > 1.9) {
          println("yikes! bond length " + d + " at " + snap.time)
        }
      }
    }
  }

  // make sure none of the monomers wrapped around the system
  def testWrap(snaps: Seq[Snapshot], mols: Seq[Seq[Int]]) {
    for (s <- snaps;
       i <- mols.flatten;
         if s.ix(i) != 0 || s.iy(i) != 0 || s.iz(i) != 0) {
      println("yikes! i="+i+", "+s.ix(i)+" "+s.iy(i)+" "+s.iz(i))
      System.exit(1)
    }
  }

  def centerOfMass(s: Snapshot, elems: Seq[Int]) = {
    var xAcc = 0.
    var yAcc = 0.
    var zAcc = 0.
    for (i <- elems) {
      xAcc += s.x(i)
      yAcc += s.y(i)
      zAcc += s.z(i)
    }
    val n = elems.length
    Vec3(xAcc/n, yAcc/n, zAcc/n)
  }

  def radiusOfGyrationSquared(s: Snapshot, elems: Seq[Int]) = {
    val cm = centerOfMass(s, elems)

    var x2Acc = 0.
    var y2Acc = 0.
    var z2Acc = 0.
    for (i <- elems) {
      x2Acc += sqr(s.x(i)-cm.x) 
      y2Acc += sqr(s.y(i)-cm.y) 
      z2Acc += sqr(s.z(i)-cm.z)
    }
    val n = elems.length
    (x2Acc/n, y2Acc/n, z2Acc/n)    
  }

  // calculate radius of gyration for entire brush
  def entireRg2(snaps: Seq[Snapshot], mols: Seq[Seq[Int]]): Seq[Double] = {
    for (s <- snaps) yield {
      val (rgx2, rgy2, rgz2) = radiusOfGyrationSquared(s, mols.flatten)
      rgx2 + rgy2 + rgz2
    }
  }

  // calculate radius of gyration for individual chains
  def chainRg2(snaps: Seq[Snapshot], mols: Seq[Seq[Int]]): Seq[Double] = {
    for (s <- snaps) yield {
      mean (mols.map {m =>
        val (rgx2, rgy2, rgz2) = radiusOfGyrationSquared(s, m)
        rgx2 + rgy2 + rgz2
      })
    }
  }

  // calculate mean end to end distance squared of chains
  def endToEnd2(snaps: Seq[Snapshot], mols: Seq[Seq[Int]]): Seq[Double] = {
    for (s <- snaps) yield {
      mean (mols.map {m =>
        val (i1, i2) = (m.head, m.last)
        val x2 = sqr(s.x(i1)-s.x(i2))
        val y2 = sqr(s.y(i1)-s.y(i2))
        val z2 = sqr(s.z(i1)-s.z(i2))
        x2 + y2 + z2
      })
    }
  }

  // jusufi's definition of brush length
  def brushLen2(snaps: Seq[Snapshot], coreRadius: Double, mols: Seq[Seq[Int]]): Seq[Double] = {
    for (s <- snaps) yield {
      // cm for entire brush
      val brushCm = centerOfMass(s, mols.flatten)

      // calculate mean center to end distance squared
      val rce2 = mean(mols.map {m =>
      brushCm distance2 (s.getPoint(m.last))
      })

      // "brush length" (aka L) squared
      sqr(sqrt(rce2) - coreRadius)
    }
  }
  
  
  def tubeDensity(snaps: Seq[Snapshot], mols: Seq[Seq[Int]]) {
    val monomers = mols.flatten
    val binWidth = 0.05
    val maxRadius = sqrt(snaps(0).maxDistance2)
    val binCnt = math.ceil(maxRadius/binWidth).intValue
    
    // centers of each bin
    val positions = (0 until binCnt) map (i => (i + 0.5) * binWidth)
    
    // density for atom types (3, 4, 5), which are
    // (counterions, salt counterions, salt coions) respectively
    val densities = (0 until 3) map (_ => new Array[Double](binCnt))
    
    // fill density arrays
    for ((s,i) <- snaps.zipWithIndex) {
      if (i % 200 == 0)
        println("Processing snapshot "+i)
      
      for (ion <- 0 until s.natoms) {
        var minDist2 = java.lang.Double.MAX_VALUE
        
        for (mon <- monomers) {
          val dist2 = s.distance2(ion, mon)
          minDist2 = math.min(minDist2, dist2)
        }

        val idx = (sqrt(minDist2) / binWidth).intValue
        s.typ(ion) match {
          case 3 => densities(0)(idx) += 1
          case 4 => densities(1)(idx) += 1
          case 5 => densities(2)(idx) += 1
          case _ => ()
        }
      }
    }
    
    // normalize density arrays 
    for (d <- densities; i <- 0 until d.size)
      d(i) /= binWidth*snaps.size
    
    // build cumulative array
    val cumulatives = for (d <- densities) yield {
      val c = new Array[Double](binCnt)
      var x = 0.0
      for (i <- 0 until binCnt) {
        x += d(i) * binWidth
        c(i) = x
      }
      c
    }
    
    // write all to file
    val formatted = formatDataInColumns(
      "radius" -> positions.toArray,
      "type3" ->  densities(0),
      "type4" ->  densities(1),
      "type5" ->  densities(2),
      "type3-cum" -> cumulatives(0),
      "type4-cum" -> cumulatives(1),
      "type5-cum" -> cumulatives(2)
    )
    
    writeStringToFile(formatted, "tube.dat")
  }
  
  //
  // Fraction of 1- ions that are associated onto Z+ ions
  //
  def associatedIons(snaps: Seq[Snapshot]): Seq[Double] = {
    // mapping from atom label to atom types index
    val monomerTypes = Set(1.0, 2.0) // graftedMonomer, monomer
    val counterionType = 3.0
    val saltCounterionType = 4.0
    val saltCoionType = 5.0
    
    val types = snaps(0).typ
    val monomers        = types.indices.filter { monomerTypes contains types(_) }
    val counterions     = types.indices.filter { counterionType == types(_) }
    val saltCounterions = types.indices.filter { saltCounterionType == types(_) }
    val saltCoions      = types.indices.filter { saltCoionType == types(_) }
    
    for ((s,i) <- snaps.zipWithIndex) yield {
      if (i % 200 == 0)
        println("Processing snapshot "+i)
      
      // is distance of atom(i) to any atom(j) less than threshold?
      def withinDistance(i: Int, js: Seq[Int], threshold: Double) =
        js exists { s.distance2(i, _) < threshold*threshold }
      
      def minDistance(is: Seq[Int], js: Seq[Int]) =
        (for (i <- is; j <- js) yield s.distance2(i,j)).min
      
      // find saltCounterions that are not too close to brush
      val freeSaltCounterions = saltCounterions filterNot { withinDistance(_, monomers, 2) }
      val boundSaltCounterions = saltCounterions filterNot (freeSaltCounterions.contains)
      
      // return fraction of saltCoions that are 'attached' to a free saltCounterion
      val associatedIons = saltCoions filter { withinDistance(_, freeSaltCounterions, 1.2) }
      
//      associatedIons.size.toDouble / saltCoions.size
      freeSaltCounterions.size.toDouble / saltCounterions.size
    }
  }

  // a note about error propagation
  // (from http://en.wikipedia.org/wiki/Propagation_of_uncertainty)
  //
  // given (x = a +- e) and (y = c x) => (y = c a +- c e)
  // given (x = a +- e) and (y = x^c) => (y = a^c +- c e (a^c / a)
  //

  def go(chainLength: Int, numChains: Int, coreRadius: Double) {
    val mols = molecules(chainLength, numChains)
    val snaps = LammpsParser.readLammpsDump("dump.dat") filter {_.time > 2000000}
    snaps.foreach {s => fixPolymerContinuity(s, mols)}
    LammpsParser.weaveThermoData(snaps, LammpsParser.readLammpsThermo("log.lammps"))
    
    // todo: return string and write to file here
    tubeDensity(snaps, mols)
    
    val formatted = formatDataInColumns(
      ("time", snaps.map(_.time.toDouble).toArray),
      ("pressure", snaps.map(_.thermo.pressure).toArray),
      ("Rg^2", entireRg2(snaps, mols).toArray),
      ("cRg^2", chainRg2(snaps, mols).toArray),
      ("Re^2", endToEnd2(snaps, mols).toArray),
      ("L^2", brushLen2(snaps, coreRadius, mols).toArray),
      ("AssociatedIons", associatedIons(snaps).toArray)
    )
    writeStringToFile(formatted, "results.dat")
//    print(formatted)

    //  val elapsedTime = timeArray(timeArray.length-1) - timeArray(0)
    //  val fftData = fft1d_continuous(rg2Array, elapsedTime).copyData().columns()
    //  writeStringToFile(formatDataInColumns("# k |FT[Rg]|^2", fftData.toSeq), "gyrfft.dat")
  }
}
