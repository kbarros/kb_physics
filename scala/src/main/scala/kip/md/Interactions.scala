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

  /** Returns force on atom 'a' due to atom 'b'.  The force on atom 'b' is
   *  the negative return value. */
  def force(that: T, a: Atom, b: Atom): Vec3

  def potential(that: T, a: Atom, b: Atom): Double
  
  def cutoff(that: T): Double
}

trait Interaction3 {
  def potential(a: Atom, b: Atom, c: Atom): Double
  def force(a: Atom, b: Atom, c: Atom): (Vec3, Vec3, Vec3)
  def cutoff(that1: Interaction3, that2: Interaction3): Double
}


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
                   val sigma_cutoff:Double=3.0) extends Interaction2 {
  type T = LennardJones
  
  def sigma_effective(that: LennardJones) = (sigma + that.sigma)/2


  override def compatibleInteractions(is: Traversable[Interaction2]): Traversable[LennardJones] = {
    is.collect { case i: LennardJones => i }
  }
  
  override def potential(that: LennardJones, a: Atom, b: Atom): Double = {
    val r2 = sqr(a.x-b.x) + sqr(a.y-b.y) + sqr(a.z-b.z)
    val sig = sigma_effective(that)
    if (r2 > sqr(sig*LennardJones.cutoff_sigma))
      0
    else {
      val a2 = sqr(sig)/r2;
      val a6 = a2*a2*a2
      val a12 = a6*a6
      4*eps*((a12 - a6) - LennardJones.shift)
    }
  }
  
  override def force(that: LennardJones, a: Atom, b: Atom): Vec3 = {
    val r2 = sqr(a.x-b.x) + sqr(a.y-b.y) + sqr(a.z-b.z)
    val sig = sigma_effective(that)
    if (r2 > sqr(sig*LennardJones.cutoff_sigma))
      Vec3(0,0,0)
    else {
      val a2 = sqr(sig)/r2;
      val a6 = a2*a2*a2
      val a12 = a6*a6
      Vec3(b.x-a.x, b.y-a.y, b.z-a.z)*(24*eps*(a12 - 2*a6)/r2)
    }
  }

  override def cutoff(that: LennardJones): Double = {
    sigma_effective(that) * sigma_cutoff
  }
}

