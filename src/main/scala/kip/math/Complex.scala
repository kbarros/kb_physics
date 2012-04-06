package kip.math


object Complex {
  implicit def scalarToComplex[T <% Double](re: T) = Complex(re, 0)
  val I = Complex(0, 1)
}

case class Complex(re: Double, im: Double) {
  def +(that: Double): Complex = Complex(re+that, im)
  def -(that: Double): Complex = Complex(re-that, im)
  def *(that: Double): Complex = Complex(re*that, im*that)
  def /(that: Double): Complex = Complex(re/that, im/that)

  def +(that: Complex): Complex = Complex(re+that.re, im+that.im)
  def -(that: Complex): Complex = Complex(re-that.re, im-that.im)
  def *(that: Complex): Complex = Complex(re*that.re - im*that.im, re*that.im + im*that.re)
  def /(that: Complex): Complex = this*that.conj / that.abs2
  
  def unary_- = Complex(-re, -im)
  def abs = math.sqrt(abs2)
  def abs2 = re*re + im*im
  def conj = Complex(re, -im)
  
  override def equals(that : Any) = that match {
    case that : Complex => re == that.re && im == that.im
    case re : Double => re == re && im == 0
    case re : Int => re == re && im == 0
    case re : Short => re == re && im == 0
    case re : Long => re == re && im == 0
    case re : Float => re == re && im == 0
    case _ => false
  }
  
  override def toString = {
    val Tol = 1e-12
    if (math.abs(im) < Tol)
      re.toString
    else if (math.abs(re) < Tol)
      im+"i"
    else {
      if (im >= 0)
        re+" + "+im+"i"
      else
        re+" - "+(-im)+"i"
    }
  }
}
