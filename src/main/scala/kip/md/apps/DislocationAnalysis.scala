package kip.md.apps

import kip.util.RangeArray

class DislocationAnalysis(tmin: Double, tmax: Double, dt: Double, rmin: Double, rmax: Double, dr: Double, volume: Double) {
  val gs: RangeArray[PairCorrelation] = RangeArray.fill(tmin, tmax, dt)(new PairCorrelation(rmin, rmax, dr, volume))
  
  def accum[A](t: Double, ids1: Seq[A], ids2: Seq[A], dist: (A, A)=>Double) {
    if (gs.isDefinedAt(t))
      gs(t).accum(ids1, ids2, dist)
  }
  
  def results: Seq[(Double, RangeArray[Double])] = for (i <- gs.indices) yield (gs.elemCenterForIndex(i), gs.elems(i).normalized)
}
