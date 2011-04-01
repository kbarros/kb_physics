package kip.md

import scala.math._

import kip.math.Vec3
import kip.math.Math._


trait Interaction1 {
  def potential(world: World, a: Atom): Double
  def force(world: World, a: Atom): Vec3
}

trait Interaction2 {
  type T <: Interaction2
  
  def compatibleInteractions(is: Traversable[Interaction2]): Traversable[T]
  
  /** Returns potential for atoms (a,b) using interactions (this, bint) */
  def potential(world: World, a: Atom, bint: T, b: Atom): Double
  
  /** Returns forces on atoms (a, b) using interactions (this, bint) */
  def force(world: World, a: Atom, bint: T, b: Atom): (Vec3, Vec3)
  
  def cutoff(that: T): Double
}

trait Interaction3 {
  type T <: Interaction3

  def compatibleInteractions(is: Traversable[Interaction3]): Traversable[T]
  
  def potential(world: World, a: Atom, bint: T, b: Atom, cint:T, c: Atom): Double
  
  def force(world: World, a: Atom, bint: T, b: Atom, cint: T, c: Atom): (Vec3, Vec3, Vec3)
  
  def cutoff(that1: T, that2: T): Double
}


// If potential and forces are symmetric between objects, then it is only necessary
// to calculate each once per pair
trait PairInteraction extends Interaction2 {
  override def potential(world: World, a: Atom, bint: T, b: Atom): Double = {
    if (a.idx == b.idx) {
      println("Error: Atoms %s and %s have the same index".format(a, b))
    }
    
    if (a.idx < b.idx) {
      val r2 = world.volume.distance2(a, b)
      if (r2 < sqr(cutoff(bint))) pairPotential(world, a, bint, b) else 0
    }
    else 0
  }
  
  override def force(world: World, a: Atom, bint: T, b: Atom): (Vec3, Vec3) = {
    if (a.idx == b.idx) {
      println("Error: Atoms %s and %s have the same index".format(a, b))
    }
    
    if (a.idx < b.idx) {
      val r2 = world.volume.distance2(a, b)
      if (r2 < sqr(cutoff(bint))) {
        pairForce(world, a, bint, b) 
      }
      else (Vec3.zero, Vec3.zero)
    }
    else (Vec3.zero, Vec3.zero)
  }


  def pairPotential(world: World, a: Atom, bint: T, b: Atom): Double
  
  def pairForce(world: World, a: Atom, bint: T, b: Atom): (Vec3, Vec3)
}


// #### LENNARD JONES #####


class LennardJones(val eps:Double=1.0,
                   val sigma_local:Double=1.0,
                   val scaled_cutoff:Double=3.0) extends PairInteraction {
  type T = LennardJones
  
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
  
  override def pairPotential(world: World, a: Atom, bint: T, b: Atom): Double = {
    val r2 = world.volume.distance2(a, b)
    val a2 = sqr(sigma(bint))/r2;
    val a6 = a2*a2*a2
    val a12 = a6*a6
    4*eps*((a12 - a6) - shift)
  }
  
  override def pairForce(world: World, a: Atom, bint: T, b: Atom): (Vec3, Vec3) = {
    val d = world.volume.displacement(a, b)
    val r2 = d.norm2
    val a2 = sqr(sigma(bint))/r2
    val a6 = a2*a2*a2
    val a12 = a6*a6
    val f = 24*eps*(2*a12 - a6)/r2
    (d*(-f), d*(f))
  }

  override def cutoff(that: LennardJones): Double = {
    sigma(that) * scaled_cutoff
  }
}


// #### SOFT #####


class PairSoft(val eps:Double=1.0,
               val sigma_local:Double=1.0) extends PairInteraction {
  type T = PairSoft
  
  def sigma(that: PairSoft) = (sigma_local + that.sigma_local)/2
  
  override def compatibleInteractions(is: Traversable[Interaction2]): Traversable[PairSoft] = {
    is.collect { case i: PairSoft => i }
  }
  
  override def pairPotential(world: World, a: Atom, bint: T, b: Atom): Double = {
    val r = sqrt(world.volume.distance2(a, b))
    val alpha = Pi/sigma(bint)
    cos(alpha*r)+1
  }
  
  override def pairForce(world: World, a: Atom, bint: T, b: Atom): (Vec3, Vec3) = {
    val d = world.volume.displacement(a, b)
    val r = d.norm
    val alpha = Pi/sigma(bint)
    val f = alpha*sin(alpha*r)
    (d*(-f/r), d*(f/r))
  }

  override def cutoff(that: PairSoft): Double = {
    sigma(that)
  }
}

