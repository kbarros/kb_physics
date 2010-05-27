package kip.util

case class Vec3(x: Double, y: Double, z: Double) {
  import Util.sqr
  
  def +(that: Vec3) = Vec3(x+that.x, y+that.y, z+that.z)
  def -(that: Vec3) = Vec3(x-that.x, y-that.y, z-that.z)
  def *(a: Double) = Vec3(x*a, y*a, z*a)
  def /(a: Double) = Vec3(x/a, y/a, z/a)
  def dot(that: Vec3) = x*that.x + y*that.y + z*that.z
  def distance2(that: Vec3) = sqr(x-that.x) + sqr(y-that.y) + sqr(z-that.z)
  
  lazy val norm = math.sqrt(x*x + y*y + z*z)
  lazy val normalize = this / norm
}
