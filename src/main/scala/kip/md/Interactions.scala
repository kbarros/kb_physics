package kip.md

import scala.math._

import kip.math.mutable
import kip.math.Vec3
import kip.math.Math._


trait Interaction1 {
  def potential(world: World, a: Atom): Double
  def accumForce(world: World, a: Atom, accum: mutable.Vec3)
}

trait Interaction2 {
  type T <: Interaction2
  
  val symmetric: Boolean
  
  def compatibleInteractions(is: Traversable[Interaction2]): Traversable[T]
  
  /** Potential energy between atoms a and b */
  def potential(world: World, a: Atom, bint: T, b: Atom): Double
  
  /** Force on atom a given interaction bint with atom b */
  def accumForce(world: World, a: Atom, bint: T, b: Atom, accum: mutable.Vec3)
  
  def cutoff(that: T): Double
}

trait Interaction3 {
  type T <: Interaction3

  def compatibleInteractions(is: Traversable[Interaction3]): Traversable[T]
  
  def potential(world: World, a: Atom, bint: T, b: Atom, cint:T, c: Atom): Double
  
  def accumForce(world: World, a: Atom, bint: T, b: Atom, cint: T, c: Atom, accum: mutable.Vec3)
  
  def cutoff(that1: T, that2: T): Double
}


// ##### LJ WALL #####

class LJWall(val pos: Vec3,
             val normal: Vec3,
             val eps: Double = 1.0,
             val sigma: Double = 1.0,
             val scaled_cutoff: Double = pow(2, 1.0/6)) extends Interaction1 {
  var netForceX = 0d
  var netForceY = 0d
  var netForceZ = 0d
  
  val shift = {
    val a2 = 1/sqr(scaled_cutoff)
    val a6 = a2*a2*a2
    val a12 = a6*a6
    a12 - a6
  }
  
  def potential(world: World, a: Atom): Double = {
    val delta = (pos - a.pos).projectOnto(normal)
    val r2 = delta.norm2
    if (r2 < sqr(cutoff)) {
      val a2 = sqr(sigma)/r2;
      val a6 = a2*a2*a2
      val a12 = a6*a6
      4*eps*((a12 - a6) - shift)
    } else {
      0
    }
  }
  
  def accumForce(world: World, a: Atom, accum: mutable.Vec3) {
    val delta = (pos - a.pos).projectOnto(normal)
    val r2 = delta.norm2
    if (r2 < sqr(cutoff)) {
      val a2 = sqr(sigma)/r2
      val a6 = a2*a2*a2
      val a12 = a6*a6
      val f = 24*eps*(2*a12 - a6)/r2
      accum -= delta*f
    }
  }

  def cutoff: Double = {
    sigma * scaled_cutoff
  }
}


// ##### LENNARD JONES #####


class LennardJones(val eps:Double=1.0,
                   val sigma_local:Double=1.0,
                   val scaled_cutoff:Double=3.0) extends Interaction2 {
  type T = LennardJones
  val symmetric = true
  
  val shift = {
    val a2 = 1/sqr(scaled_cutoff)
    val a6 = a2*a2*a2
    val a12 = a6*a6
    a12 - a6
  }
  
  def sigma(that: LennardJones) = (sigma_local + that.sigma_local)/2
  
  override def compatibleInteractions(is: Traversable[Interaction2]): Traversable[LennardJones] = {
    is.collect { case i: LennardJones => i }
  }
  
  override def potential(world: World, a: Atom, bint: T, b: Atom): Double = {
    if (a.idx != b.idx) {
      val r2 = world.volume.distance2(a, b)
      if (r2 < sqr(cutoff(bint))) {
        val a2 = sqr(sigma(bint))/r2;
        val a6 = a2*a2*a2
        val a12 = a6*a6
        4*eps*((a12 - a6) - shift)
      }
      else 0
    }
    else 0
  }
  
  override def accumForce(world: World, a: Atom, bint: T, b: Atom, accum: mutable.Vec3) {
    if (a.idx != b.idx) {
      val dx = world.volume.deltaX(a, b)
      val dy = world.volume.deltaY(a, b)
      val dz = world.volume.deltaZ(a, b)
      val r2 = dx*dx + dy*dy + dz*dz
      if (r2 < sqr(cutoff(bint))) {
        val a2 = sqr(sigma(bint))/r2
        val a6 = a2*a2*a2
        val a12 = a6*a6
        val f = 24*eps*(2*a12 - a6)/r2
        accum.x -= f*dx
        accum.y -= f*dy
        accum.z -= f*dz
      }
    }
  }

  override def cutoff(that: LennardJones): Double = {
    sigma(that) * scaled_cutoff
  }
}


// ##### SOFT #####


class PairSoft(val eps:Double=1.0,
               val sigma_local:Double=1.0) extends Interaction2 {
  type T = PairSoft
  val symmetric = true
  
  def sigma(that: PairSoft) = (sigma_local + that.sigma_local)/2
  
  override def compatibleInteractions(is: Traversable[Interaction2]): Traversable[PairSoft] = {
    is.collect { case i: PairSoft => i }
  }
  
  override def potential(world: World, a: Atom, bint: T, b: Atom): Double = {
    if (a.idx != b.idx) {
      val r2 = world.volume.distance2(a, b)
      if (r2 < sqr(cutoff(bint))) {
        val alpha = Pi/sigma(bint)
        cos(alpha*sqrt(r2))+1
      }
      else 0
    }
    else 0
  }
  
  override def accumForce(world: World, a: Atom, bint: T, b: Atom, accum: mutable.Vec3) {
    if (a.idx != b.idx) {
      val dx = world.volume.deltaX(a, b)
      val dy = world.volume.deltaY(a, b)
      val dz = world.volume.deltaZ(a, b)
      val r2 = dx*dx + dy*dy + dz*dz
      if (r2 < sqr(cutoff(bint))) {
        val r = sqrt(dx*dx + dy*dy + dz*dz)
        val alpha = Pi/sigma(bint)
        val f = alpha*sin(alpha*r)/r
        accum.x -= f*dx
        accum.y -= f*dx
        accum.z -= f*dx
      }
    }
  }

  override def cutoff(that: PairSoft): Double = {
    sigma(that)
  }
}

