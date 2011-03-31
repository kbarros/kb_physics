package kip.md

import kip.math.Math._
import kip.math._
import kip.graphics._
import scala.math._
import scala.collection.mutable.ArrayBuffer


trait Volume {
  def bounds: Bounds3d
  def wrapAtoms(atoms: Seq[Atom])
  def distance2(a: Atom, b: Atom): Double
  def buildCells(atomsPerCell: Int, atoms: Seq[Atom])
  def atomsInRange(a: Atom, range: Double): ArrayBuffer[Atom]
}


class Cuboid(lx: Double, ly: Double, lz: Double, periodic: Boolean = false) extends Volume {
  var _grid: PointGrid2d[Atom] = _

  def bounds: Bounds3d = Bounds3d(Vec3.zero, Vec3(lx, ly, lz))

  def wrapAtoms(atoms: Seq[Atom]) {
    def wrap(x: Double, wx: Int): (Double, Int) = {
      val L = _grid.L
      if (x >= 2*L || x < -L) {
        println("Simulation exploded")
        exit(-1)
      }
      if (x >= _grid.L) (x-L, wx+1)
      else if (x < 0) (x+L, wx-1)
      else (x, wx)
    }
    
    for (a <- atoms) {
      val (x, wx) = wrap(a.x, a.wx)
      val (y, wy) = wrap(a.y, a.wy)
      val (z, wz) = wrap(a.z, a.wz)
      a.x = x
      a.y = y
      a.z = z
      a.wx = wx
      a.wy = wy
      a.wz = wz
    }
  }
  
  def distance2(a: Atom, b: Atom): Double = {
    sqr(b.x-a.x) + sqr(b.y-a.y) + sqr(b.z-a.z)
  }
  
  def buildCells(atomsPerCell: Int, atoms: Seq[Atom]) {
    val cols = max(1, sqrt(atoms.size / atomsPerCell).toInt)
    if (_grid == null || _grid.cols != cols)
      _grid = new PointGrid2d(lx, cols, periodic)
    _grid.loadPoints(atoms)
  }
  
  def atomsInRange(a: Atom, range: Double): ArrayBuffer[Atom] = _grid.pointOffsetsWithinRange(a, range)
}
