package kip.math

case class Vec3(x: Double, y: Double, z: Double) {
  import Math._
  
  def norm2: Double = x*x + y*y + z*z
  def norm:  Double = math.sqrt(norm2)
  
  def normalize: Vec3 = this / norm
  def unary_- :  Vec3 = new Vec3(-x, -y, -z)
  
  def *(a: Double): Vec3 = Vec3(x*a, y*a, z*a)
  def /(a: Double): Vec3 = Vec3(x/a, y/a, z/a)
  
  def +(v: Vec3) = Vec3(x+v.x, y+v.y, z+v.z)
  def -(v: Vec3) = Vec3(x-v.x, y-v.y, z-v.z)
  
  def dot(v: Vec3) = x*v.x + y*v.y + z*v.z
  def distance2(v: Vec3) = sqr(x-v.x) + sqr(y-v.y) + sqr(z-v.z)

  def rotateBy(q: Quaternion): Vec3 = {
    assert(false, "must check")
    val vp = q.conj * Quaternion.fromVec3(this) * q
    Vec3(vp.x, vp.y, vp.z)
  }
}
