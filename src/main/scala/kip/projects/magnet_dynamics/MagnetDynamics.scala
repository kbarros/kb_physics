package kip.projects.magnet_dynamics

import scala.util.Random
import scala.math._
import kip.math.Math._
import kip.math.Vec3


class VectorField(val N: Int) {
  val x = Array.fill(N)(0.0)
  val y = Array.fill(N)(0.0)
  val z = Array.fill(N)(1.0)
  
  def copyTo(dest: VectorField) {
    assert(N == dest.N)
	Array.copy(x, 0, dest.x, 0, N)
	Array.copy(y, 0, dest.y, 0, N)
	Array.copy(z, 0, dest.z, 0, N)
  }
  
  def norm(i: Int) = math.sqrt(x(i)*x(i) + y(i)*y(i) + z(i)*z(i))
  
  def apply(i: Int) = Vec3(x(i), y(i), z(i)) 
}


trait Lattice {
  def N: Int
  def energy(s: VectorField): Double  
  def magneticField(s: VectorField, b: VectorField)
}


class MagnetDynamics(var T: Double, var alpha: Double, var dt: Double, var lattice: Lattice) {
  val N = lattice.N
  val rand = new Random(System.currentTimeMillis())
  
  val s  = new VectorField(N)
  
  // temporary space
  val sp = new VectorField(N)
  val spp = new VectorField(N)
  val b = new VectorField(N)
  
  randomizeSpins(s)

  
  def ferromagetizeSpins(s: VectorField) {
    for (i <- 0 until N) {
      s.x(i) = 0
      s.y(i) = 0
      s.z(i) = 1
    }
    normalizeSpins(s)
  }
  
  def randomizeSpins(s: VectorField) {
    for (i <- 0 until N) {
      s.x(i) = rand.nextGaussian()
      s.y(i) = rand.nextGaussian()
      s.z(i) = rand.nextGaussian()
    }
    normalizeSpins(s)
  }
  
  def normalizeSpins(s: VectorField) {
    for (i <- 0 until N) {
      val len = s.norm(i)
      s.x(i) /= len
      s.y(i) /= len
      s.z(i) /= len
    }
  }
  
  // accumulate timestep into sp:
  //   sp += scale * [dt s x (-b[s] - alpha s x b[s]) + sqrt(dt) ... ] 
  def accumEuler(s: VectorField, sp: VectorField, scale: Double) {
    lattice.magneticField(s, b)
    
    // diffusion coefficient
    val D = (alpha / (1 - alpha*alpha)) * T
    
    for (i <- 0 until N) {      
      val S = s(i)
      val B = b(i) // = - dH / dS
      val xi = Vec3(rand.nextGaussian(), rand.nextGaussian(), rand.nextGaussian())
      
      var dS = (S cross (-B  - (S cross B)  * alpha)) * dt +
               (S cross (-xi + (S cross xi) * alpha)) * (sqrt(2*D*dt))
      
      sp.x(i) += scale * dS.x
      sp.y(i) += scale * dS.y
      sp.z(i) += scale * dS.z
    }
  }
  
  // a bad integrator, for reference
  def stepEuler() {
    s.copyTo(sp)
	accumEuler(s, sp, 1.0)
	normalizeSpins(sp)
	sp.copyTo(s)
  }
  
  // a good integrator
  def stepHeun() {
	s.copyTo(sp)
	accumEuler(s, sp, 1.0)
	s.copyTo(spp)
	accumEuler(s,  spp, 0.5)
	accumEuler(sp, spp, 0.5)
	normalizeSpins(spp)
	spp.copyTo(s)
  }
}

//      /.\*/.\*/.\*
//     o - o - o -
//    /.\*/.\*/.\*      y
//   o - o - o -     /
//  /.\*/.\*/.\*    /
// o - o - o -     ----> x
//
// layer order: o: z=0
//              .: z=1
//              *: z=2

class TriLattice(val Lx: Int, val Ly: Int, Lz: Int,
                 var J1: Double, var J2: Double, var J3: Double, var Jz: Double,
                 var Dx: Double, var Dz: Double,
                 var Hz: Double) extends Lattice {
  val N = Lx*Ly
  
  val j1_neighbors = Array(
    //         o - o - o - o - o
    //        / \ / \ / \ / \ /
    //       o -[2]-[1] - o - o
    //      / \ / \ / \ / \ /
    //     o -[3]- * - [0] - o
    //    / \ / \ / \ / \ /
    //   o - o -[4]-[5]- o
    //  / \ / \ / \ / \ /
    // o - o - o - o - o
    // 
    Array( 1, 0, 0),
    Array( 0, 1, 0),
    Array(-1, 1, 0),
    Array(-1, 0, 0),
    Array( 0,-1, 0),
    Array( 1,-1, 0)
  )

  val j2_neighbors = Array(
    //         o -[1]- o - o - o
    //        / \ / \ / \ / \ /
    //      [2]- . - . -[0]- o
    //      / \ / \ / \ / \ /
    //     o - . - * - . - o
    //    / \ / \ / \ / \ /
    //   o -[3]- . - . -[5]
    //  / \ / \ / \ / \ /
    // o - o - o -[4]- o
    // 
    Array( 1, 1, 0),
    Array(-1, 2, 0),
    Array(-2, 1, 0),
    Array(-1,-1, 0),
    Array( 1,-2, 0),
    Array( 2,-1, 0)
  )
  
  val j3_neighbors = Array(
    //        [2]- . -[1]- o - o
    //        / \ / \ / \ / \ /
    //       . - . - . - . - o
    //      / \ / \ / \ / \ /
    //    [3]- . - * - . -[0]
    //    / \ / \ / \ / \ /
    //   o - . - . - . - .
    //  / \ / \ / \ / \ /
    // o - o -[4]- . -[5]
    // 
    Array( 2, 0, 0),
    Array( 0, 2, 0),
    Array(-2, 2, 0),
    Array(-2, 0, 0),
    Array( 0,-2, 0),
    Array( 2,-2, 0)
  )
  
  val jz_neighbors = Array(
    // Increasing z increases (x,y) origin of plane
    //
    // Z = +1 (upper plane)
    //     . - . - .
    //    /1\ /0\ /   cell 0 corresponds to delta (0,0,1)
    //   . - * - .
    //  / \ /2\ /
    // . - . - .
    // 
    // Z = -1 (lower plane)
    //     . - . - .
    //    / \5/ \ /   cell 3 corresponds to delta (0,0,-1)
    //   . - * - .
    //  / \3/ \4/
    // . - . - .
    // 
    Array(  0, 0, 1),
    Array( -1, 0, 1),
    Array(  0,-1, 1),
    
    Array( 0, 0, -1),
    Array( 1, 0, -1),
    Array( 0, 1, -1)
  )
  
  def couplings = Array((J1, j1_neighbors), (J2, j2_neighbors), (J3, j3_neighbors), (Jz, jz_neighbors))
  
  def latticeCoords(i: Int) = (i%Lx, (i/Lx)%Ly, (i/(Lx*Ly)))
  
  def latticeIndexUnchecked(x: Int, y: Int, z: Int) = {
    x + Lx*y + Lx*Ly*z
  }

  def latticeIndexPeriodic(x: Int, y: Int, z: Int): Int = {
    latticeIndexUnchecked((x+Lx)%Lx, (y+Ly)%Ly, (z+Lz)%Lz)
  }
  
  def latticeIndexOpenXY(x: Int, y: Int, z: Int): Int = {
    if (x < 0 || x >= Lx || y < 0 || y >= Ly) {
      -1
    }
    else {
      latticeIndexUnchecked(x, y, (z+Lz)%Lz)
    }
  }
  
  def energy(s: VectorField) = {
    var acc = 0.0
    for (i <- 0 until N) {
      val (x, y, z) = latticeCoords(i)

      // external field
      acc += - Hz * s.z(i)
      
      // anisotropies
      acc += - Dx * s.x(i) * s.x(i)
      acc += - Dz * s.z(i) * s.z(i)
      
      for ((jmag, neigh) <- couplings) {
        for (Array(dx, dy, dz) <- neigh) {
          val j = latticeIndexOpenXY(x+dx, y+dy, z+dz)
          if (j >= 0)
            acc += - 0.5 * jmag * (s(i) dot s(j)) 
        }
      }
    }
    acc
  }
  
  def magneticField(s: VectorField, b: VectorField) {
    for (i <- 0 until N) {
      b.x(i) = 0
      b.y(i) = 0
      b.z(i) = 0
      
      val (x, y, z) = latticeCoords(i)
      
      // external field
      b.z(i) += Hz

      // anisotropies
      b.x(i) += 2 * Dx * s.x(i)
      b.z(i) += 2 * Dz * s.z(i)
      
      for ((jmag, neigh) <- couplings) {
        for (Array(dx, dy, dz) <- neigh) {
          // val j = latticeIndexPeriodic(x+dx, y+dy, z+dz)
          val j = latticeIndexOpenXY(x+dx, y+dy, z+dz)
          if (j >= 0) {
            b.x(i) += jmag * s.x(j)
            b.y(i) += jmag * s.y(j)
            b.z(i) += jmag * s.z(j)
          }
        }
      }
    }
  }
}
