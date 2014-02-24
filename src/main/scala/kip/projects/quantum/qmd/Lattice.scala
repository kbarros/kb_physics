package kip.projects.quantum.qmd

import kip.math.Vec3

trait Lattice {
  def numAtoms: Int
  def initialPositions(i: Int): Array[Vec3]
  def grouping(i: Int, x: Array[Vec3], s: Int): Int
  
  def displacement(r1: Vec3, r2: Vec3) = r2 - r1
  
  def neighbors(i: Int, x: Array[Vec3], rcut: Double): Array[Int] = {
    (for (j <- 0 until numAtoms;
         if (displacement(x(i), x(j)).norm2 < rcut*rcut)) yield j).toArray
  }
}


class LinearChain(val numAtoms: Int, spacing: Double) extends Lattice {
  def initialPositions(i: Int) = {
    Array.tabulate(numAtoms) { i => Vec3(i, 0, 0)*spacing }
  }
  
  def grouping(i: Int, x: Array[Vec3], s: Int) = {
    require(numAtoms % s == 0)
    i % s
  }
}

class SquareLattice(lx: Int, ly: Int, spacing: Double) extends Lattice {
  val numAtoms = lx*ly
  
  def latticeCoords(i: Int) = {
    (i%lx, i/lx)
  }
  
  def initialPositions(i: Int) = {
    val (x, y) = latticeCoords(i)
    Array.tabulate(numAtoms) { i => Vec3(x, y, 0)*spacing }
  }
  
//  def randomizePositions(spacing: Double, disorder: Double) {
//    for (i <- 0 until numAtoms) {
//      val (x, y) = latticeCoords(i)
//      val dx = disorder*rand.nextGaussian()
//      val dy = if (ly == 1) 0 else disorder*rand.nextGaussian()
//      pos(i) = Vec3(x+0.5*y+dx, y*sqrt(3.0)/2.0+dy, 0)*spacing
//    }
//  }
  
  def grouping(i: Int, x: Array[Vec3], s: Int): Int = {
    val len = math.sqrt(s).toInt
    require(len*len == s)
    val (x, y) = latticeCoords(i)
    val xp = x % len
    val yp = y % len
    xp + yp*len
  }

}