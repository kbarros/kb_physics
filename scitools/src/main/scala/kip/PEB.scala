package kip

import Util._
import java.lang.Math._


object PEB {
    
    def molecules(chainLength: Int, numChains: Int): Seq[Seq[Int]] =
	for (i <- 0 until numChains)
	    yield (i*chainLength) until ((i+1)*chainLength)
    
    
    // make sure none of the monomers wrapped around the system
    def testWrap(snaps: Seq[Snapshot], mols: Seq[Seq[Int]]) {
	for (s <- snaps;
	     i <- mols.toList.flatten[Int]
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
	    val (rgx2, rgy2, rgz2) = radiusOfGyrationSquared(s, mols.toList.flatten[Int])
	    rgx2 + rgy2 + rgz2
	}
    }
    
    // calculate mean end to end distance squared of chains
    def endToEnd2(snaps: Seq[Snapshot], mols: Seq[Seq[Int]]): Seq[Double] = {
	for (s <- snaps) yield {
	    average (mols.map {m =>
		val (i1, i2) = (m.first, m.last)
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
	    val brushCm = centerOfMass(s, mols.toList.flatten[Int]) 
	    
	    // calculate mean center to end distance squared
	    val rce2 = average(mols.map {m =>
		brushCm distance2 (s.getPoint(m.last))
	    })
	    
	    // "brush length" (aka L) squared
	    sqr(sqrt(rce2) - coreRadius)
	}
    }

    def go(chainLength: Int, numChains: Int, coreRadius: Double) {
	val mols = molecules(chainLength, numChains)
	val snaps = LammpsParser.readLammpsDumpPartial("dump.dat", 2000)

	val times = snaps.map(_.time)
	val rg2s = entireRg2(snaps, mols)
	val ee2s = endToEnd2(snaps, mols)
	val brushLen2s = PEB.brushLen2(snaps, coreRadius, mols)

	val formatted = formatDataInColumns(("time", times.toArray),
					    ("Rg^2", rg2s.toArray),
					    ("L^2", brushLen2s.toArray),
					    ("Re^2", ee2s.toArray))
	writeStringToFile(formatted, "gyr.dat")

//	val elapsedTime = timeArray(timeArray.length-1) - timeArray(0)
//	val fftData = fft1d_continuous(rg2Array, elapsedTime).copyData().columns()
//	writeStringToFile(formatDataInColumns("# k |FT[Rg]|^2", fftData.toSeq), "gyrfft.dat")
    }
}
