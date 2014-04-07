package kip.projects.quantum.qmd

import kip.math.Vec3

object Lattice {
  def periodicWrap(dx: Double, L: Double) = {
    val dxp = dx % L
    if (dxp > L/2)
      dxp-L
    else if (dxp < -L/2)
      dxp+L
    else
      dxp
  }
}


trait Lattice {
  def numAtoms: Int
  def boundsLow: Vec3
  def boundsHigh: Vec3
  def displacement(r1: Vec3, r2: Vec3): Vec3
  def initialPositions: Array[Vec3]
  def grouping(i: Int, x: Array[Vec3], s: Int): Int
  
  def neighbors(i: Int, x: Array[Vec3], rcut: Double): Array[Int] = {
    (for (j <- 0 until numAtoms;
         if (displacement(x(i), x(j)).norm2 < rcut*rcut)) yield j).toArray
  }
}


class LinearChain(val numAtoms: Int, r0: Double, periodic: Boolean) extends Lattice {
  if (periodic) require(numAtoms > 2)
  
  val boundsLow  = Vec3(-r0, -r0, -r0)
  val boundsHigh = Vec3(numAtoms*r0, +r0, +r0)
  
  def displacement(r1: Vec3, r2: Vec3): Vec3 = {
    val delta = r2 - r1 
    if (!periodic) delta
    else {
      val dx = Lattice.periodicWrap(delta.x, numAtoms*r0)
      Vec3(dx, delta.y, delta.z)
    }
  }
  
  def initialPositions() = {
    Array.tabulate(numAtoms) { i => Vec3(i, 0, 0)*r0 }
  }
  
  def grouping(i: Int, x: Array[Vec3], s: Int) = {
    require(numAtoms % s == 0)
    i % s
  }
}

class SquareLattice(lx: Int, ly: Int, r0: Double, periodic: Boolean) extends Lattice {
  if (periodic) require(lx > 2 && ly > 2)
  
  val numAtoms = lx*ly
  val boundsLow  = Vec3(-r0, -r0, -r0)
  val boundsHigh = Vec3(lx*r0, ly*r0, +r0)
  
  def latticeCoords(i: Int) = {
    (i%lx, i/lx)
  }
  
  def displacement(r1: Vec3, r2: Vec3): Vec3 = {
    val delta = r2 - r1 
    if (!periodic) delta
    else {
      val dx = Lattice.periodicWrap(delta.x, lx*r0)
      val dy = Lattice.periodicWrap(delta.y, ly*r0)
      Vec3(dx, dy, delta.z)
    }
  }
  
  def initialPositions() = {
    Array.tabulate(numAtoms) { i =>
      val (x, y) = latticeCoords(i)
      Vec3(x, y, 0)*r0
    }
  }
  
//  def randomizePositions(r0: Double, disorder: Double) {
//    for (i <- 0 until numAtoms) {
//      val (x, y) = latticeCoords(i)
//      val dx = disorder*rand.nextGaussian()
//      val dy = if (ly == 1) 0 else disorder*rand.nextGaussian()
//      pos(i) = Vec3(x+0.5*y+dx, y*sqrt(3.0)/2.0+dy, 0)*r0
//    }
//  }
  
  def grouping(i: Int, x: Array[Vec3], s: Int): Int = {
    val len = math.sqrt(s).round.toInt
    require(len*len == s)
    val (x, y) = latticeCoords(i)
    val xp = x % len
    val yp = y % len
    xp + yp*len
  }
}

class DiamondLattice(lx: Int, ly: Int, lz: Int, r0: Double, periodic: Boolean) extends Lattice {
  if (periodic) require(lx > 1 && ly > 1 && lz > 1)
  
  val atomsPerCell = 8
  val numAtoms = atomsPerCell*lx*ly*lz
  val a = (4.0 / math.sqrt(3.0)) * r0 // lattice constant
  
  val boundsLow  = Vec3(-a, -a, -a)
  val boundsHigh = Vec3(lx*r0, ly*r0, lz*r0)
  
  
  def displacement(r1: Vec3, r2: Vec3): Vec3 = {
    val delta = r2 - r1 
    if (!periodic) delta
    else {
      val dx = Lattice.periodicWrap(delta.x, lx*a)
      val dy = Lattice.periodicWrap(delta.y, ly*a)
      val dz = Lattice.periodicWrap(delta.y, lz*a)
      Vec3(dx, dy, dz)
    }
  }
  
  def latticeCoords(i: Int) = {
    val ip = i/atomsPerCell
    val x = ip % lx
    val y = (ip / lx) % ly
    val z = ip / (lx * ly)
    (x, y, z)
  }
  
  def initialPositions() = {
    Array.tabulate(numAtoms) { i =>
      val (x, y, z) = latticeCoords(i)
      val offset = Vec3(a, a, a) / 4.0
      (i % atomsPerCell) match {
        case 0 => Vec3(x, y, z)*a
        case 1 => Vec3(x+0.5, y+0.5, z)*a
        case 2 => Vec3(x+0.5, y, z+0.5)*a
        case 3 => Vec3(x, y+0.5, z+0.5)*a
        case 4 => Vec3(x, y, z)*a         + offset
        case 5 => Vec3(x+0.5, y+0.5, z)*a + offset
        case 6 => Vec3(x+0.5, y, z+0.5)*a + offset
        case 7 => Vec3(x, y+0.5, z+0.5)*a + offset
      }
    }
  }
  
  def grouping(i: Int, x: Array[Vec3], s: Int): Int = {
    val len = math.cbrt(s/atomsPerCell).round.toInt
    require(len*len*len*atomsPerCell == s)
    val (x, y, z) = latticeCoords(i)
    val xp = x % len
    val yp = y % len
    val zp = z % len
    (xp + yp*len + zp*len*len)*atomsPerCell + (i%atomsPerCell)
  }
}
