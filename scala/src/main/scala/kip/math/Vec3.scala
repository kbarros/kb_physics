package kip.math


object Vec3 {
  val zero: Vec3 = Vec3(0,0,0)
  def apply(x: Double, y: Double, z: Double): Vec3 = immutable.Vec3(x,y,z)
}


trait Vec3 {
  import Math._

  def x: Double
  def y: Double
  def z: Double
  
  def norm2: Double = x*x + y*y + z*z
  def norm:  Double = math.sqrt(norm2)
  
  def normalize: Vec3 = this / norm
  def unary_- :  Vec3 = Vec3(-x, -y, -z)
  
  def *(a: Double): Vec3 = Vec3(x*a, y*a, z*a)
  def /(a: Double): Vec3 = Vec3(x/a, y/a, z/a)
  
  def +(v: Vec3) = Vec3(x+v.x, y+v.y, z+v.z)
  def -(v: Vec3) = Vec3(x-v.x, y-v.y, z-v.z)
  def cross(v: Vec3) = Vec3(y*v.z-z*v.y, z*v.x-x*v.z, x*v.y-y*v.x)
  // def ⨯(v: Vec3) = cross(v)
  
  def dot(v: Vec3) = x*v.x + y*v.y + z*v.z
  def distance2(v: Vec3) = sqr(x-v.x) + sqr(y-v.y) + sqr(z-v.z)
  def distance(v: Vec3) = math.sqrt(distance2(v))
  
  def projectOnto(v: Vec3): Vec3 = {
    v * ((v dot this) / v.norm2)
  }
  
  def rotateBy(q: Quaternion): Vec3 = {
    assert(false, "must check")
    val vp = q.conj * Quaternion.fromVec3(this) * q
    Vec3(vp.x, vp.y, vp.z)
  }
  
  def toMutable = mutable.Vec3(x, y, z)
}
