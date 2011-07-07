package kip.math
package linalg


object Scalar {
  implicit object RealDblTC extends RealDblTC
  implicit object RealFltTC extends RealFltTC
  implicit object ComplexTC extends ComplexTC
}

trait Scalar[@specialized(Float, Double) T] {
  def add(a: T, b: T): T
  def sub(a: T, b: T): T
  def mul(a: T, b: T): T
  def div(a: T, b: T): T
  def negate(a: T): T
  def zero: T
}

trait RealDblTC extends Scalar[Double] {
  def add(a: Double, b: Double): Double = a + b
  def sub(a: Double, b: Double): Double = a - b
  def mul(a: Double, b: Double): Double = a * b
  def div(a: Double, b: Double): Double = a / b
  def negate(a: Double): Double = -a
  def zero: Double = 0.0
}

trait RealFltTC extends Scalar[Float] {
  def add(a: Float, b: Float): Float = a + b
  def sub(a: Float, b: Float): Float = a - b
  def mul(a: Float, b: Float): Float = a * b
  def div(a: Float, b: Float): Float = a / b
  def negate(a: Float): Float = -a
  def zero: Float = 0.0f
}

trait ComplexTC extends Scalar[Complex] {
  def add(a: Complex, b: Complex): Complex = a + b
  def sub(a: Complex, b: Complex): Complex = a - b
  def mul(a: Complex, b: Complex): Complex = a * b
  def div(a: Complex, b: Complex): Complex = a / b
  def negate(a: Complex): Complex = -a
  def zero: Complex = 0
}

/*
 * TODO: move into Scalar
object ScalarOps {
  trait ScalarOps[@specialized(Float, Double) A] {
    val lhs: A
    val n: Scalar[A]
    def +(rhs:A) = n.add(lhs, rhs)
    def -(rhs:A) = n.sub(lhs, rhs)
    def *(rhs:A) = n.mul(lhs, rhs)
    def /(rhs:A) = n.div(lhs, rhs)
    def unary_-() = n.negate(lhs)
  }
  implicit def infixScalarOps[@specialized(Float, Double) A: Scalar](a: A): ScalarOps[A] = new ScalarOps[A] {
    val lhs = a
    val n = implicitly[Scalar[A]]
  }
}
*/
