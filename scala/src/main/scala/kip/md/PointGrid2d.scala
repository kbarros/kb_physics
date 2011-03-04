package kip.md

import kip.math.Math._
import kip.math.Vec3
import scala.math._
import scala.collection.mutable.ArrayBuffer


object PointGrid2d {
  def periodicOffset(L: Double, dx: Double) = {
    if (dx >  L/2)
      dx-L
    else if (dx < -L/2)
      dx+L
    else
      dx;
  }
}

class PointGrid2d[T <: Pt](L: Double, cols: Int, periodic: Boolean) {
  private var _dx = L / cols
  private val _cells: Array[ArrayBuffer[T]] = Array.fill(cols*cols) { new ArrayBuffer[T]() }
  private val tempArray = new ArrayBuffer[T]()
  private var points: Seq[T] = _
  
  loadPoints(Seq[T]())
  
  def loadPoints(points: Seq[T]) {
    this.points = points
    for (c <- _cells)
      c.clear()
    for (p <- points)
      _cells(pointToIndex(p.x, p.y)) += p
  }
  
  def pointOffsetsWithinRange(p: Pt, R: Double): ArrayBuffer[T] = {
    val imax = (R/_dx+1.0).toInt
    if (2*imax+1 > cols)
      return pointOffsetsWithinRangeSlow(p, R)
    
    val x = (p.x+L)%L
    val y = (p.y+L)%L
    val index = pointToIndex(x, y)
    val i1 = index%cols
    val j1 = index/cols
    
    tempArray.clear()
    val ret = tempArray
    // val ret = new ArrayBuffer[T]()
    
    val d2Cutoff = (sqr(R/_dx+sqrt(2))+1e-8).toInt
    
    for (di <- -imax to imax) {
      for (dj <- -imax to imax) {
        if (di*di + dj*dj <= d2Cutoff) {
          var i2 = i1+di
          var j2 = j1+dj
          
          if (periodic || (min(i2,j2) >= 0 && max(i2,j2) < cols)) {
            if (periodic) {
              i2 = (i2+cols)%cols
              j2 = (j2+cols)%cols
            }
	    
            // loop over explicit index i for speed
            val cell = _cells(cols*j2+i2)
            for (i <- 0 until cell.length) {
              val p2 = cell(i)
              val dx = p2.x - p.x + (i1+di-i2)*_dx
              val dy = p2.y - p.y + (j1+dj-j2)*_dx
              if (dx*dx + dy*dy < R*R)
                ret += p2
            }
          }
        }
      }
    }
    
    if (ret.size != pointOffsetsWithinRangeSlow(p, R).size)
      throw new IllegalStateException("Counting error.")
    
    ret
  }
  
  
  def pointOffsetsWithinRangeSlow(p: Pt, R: Double): ArrayBuffer[T] = {
    val ret = new ArrayBuffer[T]()
    
    for (i <- 0 until points.length) {
      var dx = p.x - points(i).x
      var dy = p.y - points(i).y
      if (periodic) {
        dx = PointGrid2d.periodicOffset(L, dx)
        dy = PointGrid2d.periodicOffset(L, dy)
      }
      if (dx*dx + dy*dy < R*R)
        ret += points(i)
    }
    ret
  }
  
  
  // rounding errors here are OK, as long as they occur in just this
  // one function
  private def pointToIndex(x: Double, y: Double): Int = {
    val i = (x/_dx).toInt
    val j = (y/_dx).toInt
    assert(i < cols && j < cols)
    return j*cols+i
  }
  
  // returns center of grid element
  private def indexToPoint(index: Int): Vec3 = {
    val i = index%cols
    val j = index/cols
    new Vec3((i+0.5)*_dx, (j+0.5)*_dx, 0)
  }
}
