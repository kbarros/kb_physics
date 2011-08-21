package kip.math


object Math {
  def sqr(x: Double) = x*x
  def hypot(x: Double, y: Double, z: Double) = math.sqrt(x*x + y*y + z*z)
  def avg(x: Double, y: Double) = 0.5*(x+y)
}
