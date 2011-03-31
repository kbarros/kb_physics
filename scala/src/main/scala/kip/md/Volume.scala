package kip.md

import kip.math.Math._
import kip.math._
import kip.graphics._
import scala.math._
import scala.collection.mutable.ArrayBuffer


trait Volume {
  def bounds: Bounds3d
  def wrapAtoms(atoms: Seq[Atom])
  def displacement(a: Atom, b: Atom): Vec3
  def distance2(a: Atom, b: Atom): Double
  def buildCells(atomsPerCell: Int, atoms: Seq[Atom])
  def atomsInRange(a: Atom, range: Double): ArrayBuffer[Atom]
}


object Volume {

  class Cuboid(lx: Double, ly: Double, lz: Double, periodic: Boolean = false) extends Volume {
    var grid: PointGrid2d[Atom] = new PointGrid2d(lx, ly, nx=1, ny=1, periodic)
    
    def bounds: Bounds3d = Bounds3d(Vec3.zero, Vec3(lx, ly, lz))

    def wrapAtoms(atoms: Seq[Atom]) {
      def wrap(l: Double, x: Double, wx: Int): (Double, Int) = {
        if (l == 0) {
          (0, 0)
        }
        else {
          if (x < -l || x >= 2*l) {
            println("Simulation exploded: coordinate %f outside range [%f, %f]".format(x, -l, 2*l))
            exit(-1)
          }
          if (x >= l) (x-l, wx+1)
          else if (x < 0) (x+l, wx-1)
          else (x, wx)
        }
      }
      
      for (a <- atoms) {
        val (x, wx) = wrap(lx, a.x, a.wx)
        val (y, wy) = wrap(ly, a.y, a.wy)
        val (z, wz) = wrap(lz, a.z, a.wz)
        a.x = x
        a.y = y
        a.z = z
        a.wx = wx
        a.wy = wy
        a.wz = wz
      }
    }
    
    def displacement(a: Atom, b: Atom): Vec3 = {
      grid.displacement(a, b)
    }

    def distance2(a: Atom, b: Atom): Double = {
      grid.distance2(a, b)
    }

    def buildCells(atomsPerCell: Int, atoms: Seq[Atom]) {
      assert(false) // need to implement 3d version
      
      // nx ny nz = atoms.size / atomsPerCell
      // nx/ny = lx/ly, etc.
      val nx = max(1, sqrt((lx/ly) * atoms.size / atomsPerCell).toInt)
      val ny = max(1, sqrt((ly/lx) * atoms.size / atomsPerCell).toInt)
      val nz = 
      if (grid.nx != nx || grid.ny != ny)
        grid = new PointGrid2d(lx, ly, nx, ny, periodic)
      grid.loadPoints(atoms)
    }
    
    def atomsInRange(a: Atom, range: Double): ArrayBuffer[Atom] = grid.pointOffsetsWithinRange(a, range)
  }
}
