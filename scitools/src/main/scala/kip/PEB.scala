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
    
    // calculate radius of gyration for entire brush
    def rg2(snaps: Seq[Snapshot], mols: Seq[Seq[Int]]): Seq[Double] = {
	for (s <- snaps) yield {
	    // calculate mean radius of gyration for polymers
	    /*
	     val rgp2 = average (mols.map {m =>	
	     val (x2, y2, z2) = s.radiusOfGyrationSquared(m)
	     x2+y2+z2
	     })
	     */
	    
	    val (rgx2, rgy2, rgz2) = s.radiusOfGyrationSquared(mols.toList.flatten[Int])
	    rgx2 + rgy2 + rgz2	    
	}
    }
    
    // jusufi's definition of brush length
    def brushLen2(snaps: Seq[Snapshot], coreRadius: Double, mols: Seq[Seq[Int]]): Seq[Double] = {
	for (s <- snaps) yield {
	    // cm for entire brush
	    val brushCm = s.centerOfMass(mols.toList.flatten[Int]) 
	    
	    // calculate mean center to end distance squared
	    val rce2 = average(mols.map {m =>
		brushCm distance2 (s.getPoint(m.last))
	    })
	    
	    // "brush length" (aka L) squared
	    sqr(sqrt(rce2) - coreRadius)
	}
    }
}
