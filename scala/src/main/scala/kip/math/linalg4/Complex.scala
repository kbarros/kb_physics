package kip.math.linalg4

case class Complexd(re: Double, im: Double) {
  def +(that: Double): Complexd = Complexd(re+that, im)
  def -(that: Double): Complexd = Complexd(re-that, im)
  def *(that: Double): Complexd = Complexd(re*that, im*that)
  def /(that: Double): Complexd = Complexd(re/that, im/that)

  def +(that: Complexd): Complexd = Complexd(re+that.re, im+that.im)
  def -(that: Complexd): Complexd = Complexd(re-that.re, im-that.im)
  def *(that: Complexd): Complexd = Complexd(re*that.re - im*that.im, re*that.im + im*that.re)
  def /(that: Complexd): Complexd = this*that.conj / that.abs2
  
  def unary_- = Complexd(-re, -im)
  def abs = math.sqrt(abs2)
  def abs2 = re*re + im*im
  def conj = Complexd(re, -im)
  
  override def equals(that : Any) = that match {
    case that : Complexd => re == that.re && im == that.im
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

case class Complexf(re: Float, im: Float) {
  def +(that: Float): Complexf = Complexf(re+that, im)
  def -(that: Float): Complexf = Complexf(re-that, im)
  def *(that: Float): Complexf = Complexf(re*that, im*that)
  def /(that: Float): Complexf = Complexf(re/that, im/that)

  def +(that: Complexf): Complexf = Complexf(re+that.re, im+that.im)
  def -(that: Complexf): Complexf = Complexf(re-that.re, im-that.im)
  def *(that: Complexf): Complexf = Complexf(re*that.re - im*that.im, re*that.im + im*that.re)
  def /(that: Complexf): Complexf = this*that.conj / that.abs2
  
  def unary_- = Complexf(-re, -im)
  def abs = math.sqrt(abs2)
  def abs2 = re*re + im*im
  def conj = Complexf(re, -im)
  
  override def equals(that : Any) = that match {
    case that : Complexf => re == that.re && im == that.im
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
