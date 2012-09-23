package kip.projects.kondo_rkky

import scala.util.Random
import scala.math._
import kip.math.fft.FFTReal
import kip.math.Math._


class Lattice3d(lx: Int, ly: Int, lz: Int) {  
  def idx(x: Int, y: Int, z: Int) = {
    mod(z, lz)*(lx*ly) + mod(y, ly)*(lx) + mod(x, lx)
  }
  
  lazy val coords: Array[(Int, Int, Int)] = {
    Array.tabulate(lx*ly*lz) { i =>
      val x = i % lx
      val y = (i / lx) % ly
      val z = i / (lx *ly)
      (x, y, z)
    }
  }
  
  lazy val periodicNeighbors: Array[Array[Int]] = {
    Array.tabulate(lx*ly*lz) { i =>
      val (x, y, z) = coords(i)
      Array(
        idx(x+1, y, z),
        idx(x-1, y, z),
        idx(x, y+1, z),
        idx(x, y-1, z),
        idx(x, y, z+1),
        idx(x, y, z-1)
      )
    }
  }
}

class DipoleSim2(val L: Int, val depth: Int, var T: Double, var H: Double, var J: Double, var anisotropy: Double, var dt: Double) {
  val Lx = L
  val Ly = L
  val Lz = L
  
  val lat = new Lattice3d(Lx, Ly, Lz)
  
  val rand = new Random(System.currentTimeMillis())
  
  val dimensions = Array(L, L, L)
  val N = L*L*L
  
  val sx = new Array[Double](N)
  val sy = new Array[Double](N)
  val sz = new Array[Double](N)
  randomizeSpins()
  
  val sbarx = new Array[Double](N)
  val sbary = new Array[Double](N)
  val sbarz = new Array[Double](N)
  
  val fft = new FFTReal(Array(Lz, Ly, Lx))

  val sqx = fft.allocFourierArray()
  val sqy = fft.allocFourierArray()
  val sqz = fft.allocFourierArray()
  
  val ks = fft.allFourierVectors
  
  // sbar_k = k (k dot s_k) / k^2
  
  def calculateSbar() {
    fft.forwardTransform(sx, sqx)
    fft.forwardTransform(sy, sqy)
    fft.forwardTransform(sz, sqz)

    for (i <- 0 until fft.nrecip/2) {
      val kz = ks(i)(0)
      val ky = ks(i)(1)
      val kx = ks(i)(2)
      
      val k2 = kx*kx + ky*ky + kz*kz + 1e-6
      val k_dot_s_re = kx * sqx(2*i+0) + ky * sqy(2*i+0) + kz * sqz(2*i+0)
      val k_dot_s_im = kx * sqx(2*i+1) + ky * sqy(2*i+1) + kz * sqz(2*i+1)
      
      sqx(2*i+0) = kx * (k_dot_s_re) / k2
      sqy(2*i+0) = ky * (k_dot_s_re) / k2
      sqz(2*i+0) = kz * (k_dot_s_re) / k2
      
      sqx(2*i+1) = kx * (k_dot_s_im) / k2
      sqy(2*i+1) = ky * (k_dot_s_im) / k2
      sqz(2*i+1) = kz * (k_dot_s_im) / k2
    }

    fft.backwardTransform(sqx, sbarx)
    fft.backwardTransform(sqy, sbary)
    fft.backwardTransform(sqz, sbarz)
    
    // nearest neighbor ferromagnetic exchange
    for (i <- 0 until Lx*Ly*Lz;
         j <- lat.periodicNeighbors(i)) {
      sbarx(i) -= J * sx(j)
      sbary(i) -= J * sy(j)
      sbarz(i) -= J * sz(j)      
    }
    
  }
  
  def randomizeSpins() {
    for (i <- 0 until N) {
      sx(i) = rand.nextGaussian()
      sy(i) = rand.nextGaussian()
      sz(i) = rand.nextGaussian()
    }
    normalizeSpins()
  }

  def normalizeSpins() {
    for (i <- 0 until L*L*L) {
      val (x, y, z) = lat.coords(i)
      if (z < depth) {
        val len = math.sqrt(sx(i)*sx(i) + sy(i)*sy(i) + sz(i)*sz(i))
        sx(i) /= len
        sy(i) /= len
        sz(i) /= len
      }
      else {
        sx(i) = 0
        sy(i) = 0
        sz(i) = 0
      }
    }
  }
  
  // sum_q chi_q |S(q)|^2
  def energy(): Double = {
    calculateSbar()
    
    var ret = 0.0
    for (i <- 0 until N) {
      ret += 0.5 * (sx(i)*sbarx(i) + sy(i)*sbary(i) + sz(i)*sbarz(i))
    }
    ret / N
  }
  
  def step() {
    calculateSbar()

    for (i <- 0 until N) {
      sx(i) -= dt * sbarx(i) + math.sqrt(dt * 2 * T) * rand.nextGaussian()
      sy(i) -= dt * sbary(i) + math.sqrt(dt * 2 * T) * rand.nextGaussian()
      sz(i) -= dt * sbarz(i) + math.sqrt(dt * 2 * T) * rand.nextGaussian()

      sz(i) += dt * (H + anisotropy * sz(i))
    }
    
    normalizeSpins()
  }
}
