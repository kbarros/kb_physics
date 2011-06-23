package kip.md

import kip.math.Math._
import kip.math.Vec3
import scala.math._
import scala.collection.mutable.ArrayBuffer


object PointGrid2d {
  trait Pt {
    def pos: Vec3
  }

  def periodicOffset(L: Double, dx: Double) = {
    if (dx > L/2)
      dx-L
    else if (dx < -L/2)
      dx+L
    else
      dx;
  }
  
  
  def cellDimensions(lx: Double, ly: Double, lz:Double, ncells:Int): (Int, Int, Int) = {
    assert (ncells > 0)
    
    def cellDimensionsAux(lx: Double, ly: Double, lz: Double): (Int, Int, Int) = {
      assert (lx >= ly && ly >= lz)
      
      def dims1d: (Int, Int, Int) = {
        (ncells, 1, 1)
      }
      
      def dims2d: (Int, Int, Int) = {
        if (ly == 0)
          dims1d
        else {
          val vol = lx*ly
          val nx = lx*sqrt(ncells/vol)
          val ny = ly*sqrt(ncells/vol)
          if (ny >= 1.5)
            (nx.round.toInt, ny.round.toInt, 1)
          else
            dims1d
        }
      }
      
      def dims3d: (Int, Int, Int) = {
        if (lz == 0)
          dims2d
        else {
          val vol = lx*ly*lz
          val nx = lx*cbrt(ncells/vol)
          val ny = ly*cbrt(ncells/vol)
          val nz = lz*cbrt(ncells/vol)
          if (nz >= 1.5)
            (nx.round.toInt, ny.round.toInt, nz.round.toInt)
          else
            dims2d
        }
      }
      
      dims3d
    }
    
    val (Seq(lxp, lyp, lzp), perm) = Seq(lx, ly, lz).zipWithIndex.sortBy(- _._1).unzip
    val (nxp, nyp, nzp) = cellDimensionsAux(lxp, lyp, lzp)
    val Seq(nx, ny, nz) = (Seq(nxp, nyp, nzp) zip perm).sortBy(_._2).unzip._1
    (nx, ny, nz)
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
      _cells(pointToIndex(p.pos.x, p.pos.y)) += p
  }
  
  def deltaX(p1: PointGrid2d.Pt, p2: PointGrid2d.Pt): Double = {
    if (periodic)
      PointGrid2d.periodicOffset(Lx, p2.pos.x - p1.pos.x)
    else
      p2.pos.x - p1.pos.x
  }

  def deltaY(p1: PointGrid2d.Pt, p2: PointGrid2d.Pt): Double = {
    if (periodic)
      PointGrid2d.periodicOffset(Ly, p2.pos.y - p1.pos.y)
    else
      p2.pos.y - p1.pos.y
  }

  def deltaZ(p1: PointGrid2d.Pt, p2: PointGrid2d.Pt): Double = {
    // if (periodic) {
    //   PointGrid2d.periodicOffset(Lz, p2.z - p1.z)
    // }
    // else
    //   p2.z - p1.z
    0
  }
  
  def distance2(p1: PointGrid2d.Pt, p2: PointGrid2d.Pt): Double = {
    sqr(deltaX(p1,p2)) + sqr(deltaY(p1,p2))
  }
  
  def pointOffsetsWithinRange(p: PointGrid2d.Pt, R: Double): ArrayBuffer[T] = {
    val imax = if (Lx == 0.0) 0 else (R/_dx+1.0).toInt
    val jmax = if (Ly == 0.0) 0 else (R/_dy+1.0).toInt
    if (2*imax+1 > nx || 2*jmax+1 > ny)
      return pointOffsetsWithinRangeSlow(p, R)
    
    val x = (p.pos.x+Lx)%Lx
    val y = (p.pos.y+Ly)%Ly
    val index = pointToIndex(x, y)
    val i1 = index%nx
    val j1 = index/nx
    
    tempArray.clear()
    val ret = tempArray
    // val ret = new ArrayBuffer[T]()

    for (di <- -imax to imax) {
      for (dj <- -jmax to jmax) {
        // two cells won't contain interacting particles if the nearest corner-to-corner
        // distance (a) is greater than (R): a > R
        // the center-to-center distance (b) minus the cell-diagonal (c) is less than the
        // corner-to-corner distance: a > b - c
        // thus, the condition, b - c > R guarantees a > R
        // for efficiency, we use the form b^2 > (R + c)^2
        val centerToCenter2 = sqr(di*_dx) + sqr(dj*_dy)
        val cellDiagonal = sqrt(_dx*_dx + _dy*_dy)
        val cellsExcluded = centerToCenter2 > sqr((R+1e-8)+cellDiagonal)
        if (!cellsExcluded) {
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
              val dx = p2.pos.x - p.pos.x + (i1+di-i2)*_dx
              val dy = p2.pos.y - p.pos.y + (j1+dj-j2)*_dy
              if (dx*dx + dy*dy < R*R)
                ret += p2
            }
          }
        }
      }
    }
    
    //if (ret.size != pointOffsetsWithinRangeSlow(p, R).size)
    //  throw new IllegalStateException("Counting error.")
    
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
    (0 <= i && i < nx) && (0 <= j && j < ny)
  }
  
  
  // rounding errors here are OK, as long as they occur in just this
  // one function
  private def pointToIndex(x: Double, y: Double): Int = {
    val i = (x/_dx).toInt
    val j = (y/_dy).toInt
    assert(validCoords(i, j))
    return j*nx+i
  }
  
  // returns center of grid element
  private def indexToPoint(index: Int): Vec3 = {
    val i = index%nx
    val j = index/ny
    Vec3((i+0.5)*_dx, (j+0.5)*_dy, 0)
  }
}
