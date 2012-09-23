package kip.projects.kondo_rkky

import scala.util.Random
import scala.math._
import kip.math.fft.FFTReal
import kip.math.Math._


class DipoleSim(val Lx: Int, val Ly: Int, val Lz: Int, var T: Double, var H: Double, var J: Double, var anisotropy: Double, var dt: Double) {
  val rand = new Random(System.currentTimeMillis())

  val N = Lx*Ly*Lz
  
  val sx = new Array[Double](N)
  val sy = new Array[Double](N)
  val sz = new Array[Double](N)
  randomizeSpins()
  
  val sbarx = new Array[Double](N)
  val sbary = new Array[Double](N)
  val sbarz = new Array[Double](N)
  
  // sbar_1 = sum_y V_{12} dot s_2
  // V = - 3 (r \tensor r) /r^5 + I / r^3
  // V dot s2 = - 3 r (r dot w2) / r^5 + s2 / r^3
  
  def calculateSbar() {
    for (x1 <- 0 until Lx;
         y1 <- 0 until Ly;
         z1 <- 0 until Lz) {
      val i1 = z1*(Lx*Ly) + y1*(Lx) + x1
      sbarx(i1) = 0
      sbary(i1) = 0
      sbarz(i1) = 0
      
      for (x2 <- 0 until Lx;
           y2 <- 0 until Ly;
           z2 <- 0 until Lz;
           i2 = z2*(Lx*Ly) + y2*(Lx) + x2;
           if (i1 != i2)) {
        val rx = x2-x1
        val ry = y2-y1
        val rz = z2-z1

        val r2 = rx*rx + ry*ry + rz*rz
        val r = math.sqrt(r2)

        val r3 = r*r*r
        val r5 = r3*r*r

        val r_dot_s2 = rx*sx(i2) + ry*sy(i2) + rz*sz(i2)

        sbarx(i1) += (- 3*rx*r_dot_s2 / r5 + sx(i2) / r3)
        sbary(i1) += (- 3*ry*r_dot_s2 / r5 + sy(i2) / r3)
        sbarz(i1) += (- 3*rz*r_dot_s2 / r5 + sz(i2) / r3)
        
        // ferromagnetic exchange interaction for nearest neighbors
        if (r2 == 1) {
          sbarx(i1) -= J*sx(i2)
          sbary(i1) -= J*sy(i2)
          sbarz(i1) -= J*sz(i2)
        }
      }
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
    for (i <- 0 until N) {
      val len = math.sqrt(sx(i)*sx(i) + sy(i)*sy(i) + sz(i)*sz(i))
      sx(i) /= len
      sy(i) /= len
      sz(i) /= len
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
