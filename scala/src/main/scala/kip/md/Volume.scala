package kip.md

import kip.math.Math._
import kip.math._
import kip.graphics._
import scala.math._
import scala.collection.mutable.ArrayBuffer


trait Volume {
  def bounds: Bounds3d
  def wrapAtoms(atoms: Seq[Atom])
  def deltaX(a: Atom, b: Atom): Double
  def deltaY(a: Atom, b: Atom): Double
  def deltaZ(a: Atom, b: Atom): Double
  def displacement(a: Atom, b: Atom): Vec3 = Vec3(deltaX(a,b), deltaY(a,b), deltaZ(a,b))
  def distance2(a: Atom, b: Atom): Double = sqr(deltaX(a,b)) + sqr(deltaY(a,b)) + sqr(deltaZ(a,b))
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
            sys.exit(-1)
          }
          if (x >= l) (x-l, wx+1)
          else if (x < 0) (x+l, wx-1)
          else (x, wx)
        }
      }
      
      for (a <- atoms) {
        val (x, wx) = wrap(lx, a.pos.x, a.wx)
        val (y, wy) = wrap(ly, a.pos.y, a.wy)
        val (z, wz) = wrap(lz, a.pos.z, a.wz)
        a.pos.x = x
        a.pos.y = y
        a.pos.z = z
        a.wx = wx
        a.wy = wy
        a.wz = wz
        
        if (!periodic && (wx != 0 || wy != 0 || wz != 0)) {
          println("Atom <%s> escaped boundary. Exiting.".format(a))
          sys.exit(-1)
        }
      }
    }
    
    def deltaX(p1: Atom, p2: Atom): Double = grid.deltaX(p1, p2)
    def deltaY(p1: Atom, p2: Atom): Double = grid.deltaY(p1, p2)
    def deltaZ(p1: Atom, p2: Atom): Double = grid.deltaZ(p1, p2)
    
    def buildCells(atomsPerCell: Int, atoms: Seq[Atom]) {
      val ncells = max(1, atoms.size / atomsPerCell)
      val (nx, ny, nz) = PointGrid2d.cellDimensions(lx, ly, lz, ncells)
      if (grid.nx != nx || grid.ny != ny)
        grid = new PointGrid2d(lx, ly, nx, ny, periodic)
      grid.loadPoints(atoms)
    }
    
    def atomsInRange(a: Atom, range: Double): ArrayBuffer[Atom] = grid.pointOffsetsWithinRange(a, range)
  }
}
