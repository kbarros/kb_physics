package kip.md

import kip.math.Math._
import kip.math.Vec3
import scala.math._
import scala.collection.mutable.ArrayBuffer


object PointGrid2d {
  trait Pt {
    var x, y: Double
    def toVec = Vec3(x, y, 0)
  }

  def periodicOffset(L: Double, dx: Double) = {
    if (dx > L/2)
      dx-L
    else if (dx < -L/2)
      dx+L
    else
      dx;
  }
}

class PointGrid2d[T <: PointGrid2d.Pt](val Lx: Double, val Ly: Double, val nx: Int, val ny: Int, val periodic: Boolean) {
  assert (nx > 0 && ny > 0)
  
  private var _dx = Lx / nx
  private var _dy = Ly / ny
  private val _cells: Array[ArrayBuffer[T]] = Array.fill(nx*ny) { new ArrayBuffer[T]() }
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

  def distance2(p1: PointGrid2d.Pt, p2: PointGrid2d.Pt): Double = {
    var dx = p2.x - p1.x
    var dy = p2.y - p1.y
    if (periodic) {
      dx = PointGrid2d.periodicOffset(Lx, dx)
      dy = PointGrid2d.periodicOffset(Ly, dy)
    }
    dx*dx + dy*dy
  }

  def displacement(p1: PointGrid2d.Pt, p2: PointGrid2d.Pt): Vec3 = {
    var dx = p2.x - p1.x
    var dy = p2.y - p1.y
    if (periodic) {
      dx = PointGrid2d.periodicOffset(Lx, dx)
      dy = PointGrid2d.periodicOffset(Ly, dy)
    }
    Vec3(dx, dy, 0)
  }
  
  def pointOffsetsWithinRange(p: PointGrid2d.Pt, R: Double): ArrayBuffer[T] = {
    val imax = if (Lx == 0.0) 0 else (R/_dx+1.0).toInt
    val jmax = if (Ly == 0.0) 0 else (R/_dy+1.0).toInt
    if (2*imax+1 > nx || 2*jmax+1 > ny)
      return pointOffsetsWithinRangeSlow(p, R)
    
    val x = (p.x+Lx)%Lx
    val y = (p.y+Ly)%Ly
    val index = pointToIndex(x, y)
    val i1 = index%nx
    val j1 = index/nx
    
    tempArray.clear()
    val ret = tempArray
    // val ret = new ArrayBuffer[T]()

    // TODO FIXME
    val d2Cutoff = (sqr(R/_dx+sqrt(2))+1e-8).toInt
    
    for (di <- -imax to imax) {
      for (dj <- -jmax to jmax) {
        if (di*di + dj*dj <= d2Cutoff) {
          var i2 = i1+di
          var j2 = j1+dj
          
          if (periodic || validCoords(i2, j2)) {
            if (periodic) {
              i2 = (i2+nx)%nx
              j2 = (j2+ny)%ny
            }
	    
            // loop over explicit index i for speed
            val cell = _cells(nx*j2+i2)
            for (i <- 0 until cell.length) {
              val p2 = cell(i)
              val dx = p2.x - p.x + (i1+di-i2)*_dx
              val dy = p2.y - p.y + (j1+dj-j2)*_dy
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
  
  
  def pointOffsetsWithinRangeSlow(p: PointGrid2d.Pt, R: Double): ArrayBuffer[T] = {
    val ret = new ArrayBuffer[T]()
    
    for (i <- 0 until points.length) {
      if (distance2(p, points(i)) < R*R)
        ret += points(i)
    }
    ret
  }
  
  
  private def validCoords(i: Int, j: Int): Boolean = {
    (0 < i && i < nx) && (0 < j && j < ny)
  }
  
  
  // rounding errors here are OK, as long as they occur in just this
  // one function
  private def pointToIndex(x: Double, y: Double): Int = {
    val i = (x/_dx).toInt
    val j = (y/_dy).toInt
    println(i + " "+j)
    assert(validCoords(i, j))
    return j*nx+i
  }
  
  // returns center of grid element
  private def indexToPoint(index: Int): Vec3 = {
    val i = index%nx
    val j = index/ny
    new Vec3((i+0.5)*_dx, (j+0.5)*_dy, 0)
  }
}
