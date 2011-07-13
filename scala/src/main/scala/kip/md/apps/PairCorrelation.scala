package kip.md.apps

import kip.md._
import kip.math.Vec3
import scala.math._
import kip.util.RangeArray


class PairCorrelation {
  // val g = RangeArray.fill(xmin=0, xmax=rmax, dx=dr)(0d)
  
  def pairCorrelation(ids1: Seq[Int], ids2: Seq[Int], dist: (Int, Int)=>Double, g: RangeArray[Double], volume: Double): Unit = {
    
    // number of unique ordered pairs between ids1 and ids2
    val uniquePairs = (ids1.size * ids2.size) - (ids1.toSet & ids2.toSet).size

    // sum over all snapshots, and all pairs of particles. distances are binned into g
    for (i1 <- ids1; i2 <- ids2; if i1 != i2) {
      val r = dist(i1, i2)
      if (g.isDefinedAt(r)) {
        // ratio of shell volume (width dr) to system volume
        val volume_fraction = 2*math.Pi*r*g.dx / volume
        // normalize pair-correlation so that it becomes unity at homogeneous density
        g(r) += 1 / (uniquePairs * volume_fraction)
      }
    }
  }
}
