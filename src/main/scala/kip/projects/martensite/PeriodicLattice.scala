package kip.projects.martensite

class PeriodicLattice(val lx: Int, val ly: Int, val lz: Int, val dx: Double) {
  val N = lx*ly*lz
  
  def coord2index(x: Int, y: Int, z: Int) = {
    val xp = (x+lx)%lx
    val yp = (y+ly)%ly
    val zp = (z+lz)%lz
    xp + yp*lx + zp*lx*ly
  }
  
  def index2coord(i: Int): (Int, Int, Int) = {
    val x = i % lx
    val y = (i % (lx*ly)) / lx
    val z = i / (lx*ly)
    (x, y, z)
  }
  
  val sE = Array.tabulate(N) { i => val (x, y, z) = index2coord(i); coord2index(x+1, y, z) }
  val sW = Array.tabulate(N) { i => val (x, y, z) = index2coord(i); coord2index(x-1, y, z) }
  val sN = Array.tabulate(N) { i => val (x, y, z) = index2coord(i); coord2index(x, y+1, z) }
  val sS = Array.tabulate(N) { i => val (x, y, z) = index2coord(i); coord2index(x, y-1, z) }
  val sU = Array.tabulate(N) { i => val (x, y, z) = index2coord(i); coord2index(x, y, z+1) }
  val sD = Array.tabulate(N) { i => val (x, y, z) = index2coord(i); coord2index(x, y, z-1) }
  
  def dX(a: Array[Double], i: Int): Double =
    (a(sE(i)) - a(sW(i))) / (2*dx)
  def dY(a: Array[Double], i: Int): Double =
    (a(sN(i)) - a(sS(i))) / (2*dx)
  def dZ(a: Array[Double], i: Int): Double =
    (a(sU(i)) - a(sD(i))) / (2*dx)
  def laplacian(a: Array[Double], i: Int): Double =
    (a(sE(i)) + a(sW(i)) + a(sN(i)) + a(sS(i)) + a(sU(i)) + a(sD(i)) - 6*a(i)) / (dx*dx)
}
