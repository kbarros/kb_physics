package kip.md

import kip.math.Vec3
import kip.math.Math._


trait Interaction1 {
  def potential(a: Atom): Double
  def force(a: Atom): Vec3
}

trait Interaction2 {
  type T <: Interaction2
  
  def compatibleInteractions(is: Traversable[Interaction2]): Traversable[T]
  
  /** Returns potential for atoms (a,b) using interactions (this, bint) */
  def potential(a: Atom, bint: T, b: Atom): Double
  
  /** Returns forces on atoms (a, b) using interactions (this, bint) */
  def force(a: Atom, bint: T, b: Atom): (Vec3, Vec3)
  
  def cutoff(that: T): Double
}

trait Interaction3 {
  type T <: Interaction3

  def compatibleInteractions(is: Traversable[Interaction3]): Traversable[T]
  
  def potential(a: Atom, bint: T, b: Atom, cint:T, c: Atom): Double
  
  def force(a: Atom, bint: T, b: Atom, cint: T, c: Atom): (Vec3, Vec3, Vec3)
  
  def cutoff(that1: T, that2: T): Double
}


// If potential and forces are symmetric between objects, then it is only necessary
// to calculate each once per pair
trait PairInteraction extends Interaction2 {

  override def potential(a: Atom, bint: T, b: Atom): Double = {
    if (a.idx < b.idx) {
      val r2 = sqr(b.x-a.x) + sqr(b.y-a.y) + sqr(b.z-a.z)
      if (r2 < sqr(cutoff(bint))) pairPotential(a, bint, b) else 0
    }
    else 0
  }
  
  override def force(a: Atom, bint: T, b: Atom): (Vec3, Vec3) = {
    if (a.idx < b.idx) {
      val r2 = sqr(b.x-a.x) + sqr(b.y-a.y) + sqr(b.z-a.z)
      if (r2 < sqr(cutoff(bint))) {
        pairForce(a, bint, b) 
      }
      else (Vec3.zero, Vec3.zero)
    }
    else (Vec3.zero, Vec3.zero)
  }


  def pairPotential(a: Atom, bint: T, b: Atom): Double
  
  def pairForce(a: Atom, bint: T, b: Atom): (Vec3, Vec3)
}


// #### LENNARD JONES #####

object LennardJones {
  val cutoff_sigma = 3 // cutoff in units of sigma
  
  val shift = {
    val a2 = 1/sqr(cutoff_sigma)
    val a6 = a2*a2*a2
    val a12 = a6*a6
    a12 - a6
  }
}


class LennardJones(val eps:Double=1.0,
                   val sigma:Double=1.0,
                   val sigma_cutoff:Double=3.0) extends PairInteraction {
  type T = LennardJones
  
  def sigma_effective(that: LennardJones) = (sigma + that.sigma)/2

  override def compatibleInteractions(is: Traversable[Interaction2]): Traversable[LennardJones] = {
    is.collect { case i: LennardJones => i }
  }
  
  override def pairPotential(a: Atom, bint: T, b: Atom): Double = {
    val r2 = sqr(b.x-a.x) + sqr(b.y-a.y) + sqr(b.z-a.z)
    val sig = sigma_effective(bint)
    val a2 = sqr(sig)/r2;
    val a6 = a2*a2*a2
    val a12 = a6*a6
    4*eps*((a12 - a6) - LennardJones.shift)
  }
  
  override def pairForce(a: Atom, bint: T, b: Atom): (Vec3, Vec3) = {
    val d = Vec3(b.x-a.x, b.y-a.y, b.z-a.z)
    val r2 = d.norm2
    val sig = sigma_effective(bint)
    val a2 = sqr(sig)/r2
    val a6 = a2*a2*a2
    val a12 = a6*a6
    val f = (24*eps*(a12 - 2*a6)/r2)
    (d*f, d*(-f))
  }

  override def cutoff(that: LennardJones): Double = {
    sigma_effective(that) * sigma_cutoff
  }
}

