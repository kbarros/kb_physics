package kip.md.apps

import kip.md._
import kip.util.RangeArray


class PairCorrelation(rmin: Double, rmax: Double, dr: Double, volume: Double) {
  var cnt: Int = 0
  val g: RangeArray[Double] = RangeArray.fill(rmin, rmax, dr)(0)
  
  def accum[A](ids1: Seq[A], ids2: Seq[A], dist: (A, A)=>Double) {
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
    cnt += 1
  }
  
  def normalized: RangeArray[Double] = g.map(x => if (x == 0) 0 else x / cnt)
}
