package kip.math


object Quaternion {
  val I = Quaternion(0, 1, 0, 0)
  val J = Quaternion(0, 0, 1, 0)
  val K = Quaternion(0, 0, 0, 1)
  
  def fromAxisAngle(x: Double, y: Double, z: Double): Quaternion = {
    val theta = math.sqrt(x*x + y*y + z*z)
    val n1 = math.cos(theta/2)
    val n2 = if (theta == 0) 0 else math.sin(theta/2)/theta
    Quaternion(n1, n2*x, n2*y, n2*z)
  }
  
  def fromVec3(v: Vec3): Quaternion = {
    Quaternion(0, v.x, v.y, v.z)
  }
}


case class Quaternion(val w: Double, val x: Double, val y: Double, val z: Double) {
  def norm2: Double = w*w + x*x + y*y + z*z
  def norm:  Double = math.sqrt(norm2)
  
  def normalize: Quaternion = this / norm
  def unary_- :  Quaternion = Quaternion(-w, -x, -y, -z)
  
  def conj:      Quaternion = Quaternion(w, -x, -y, -z)
  def inverse:   Quaternion = conj / norm2
  
  def *(a: Double): Quaternion = Quaternion(w*a, x*a, y*a, z*a)
  def /(a: Double): Quaternion = Quaternion(w/a, x/a, y/a, z/a)
  
  def *(q: Quaternion): Quaternion =
    Quaternion(
      w = + w*q.w - x*q.x - y*q.y - z*q.z,
      x = + w*q.x + x*q.w + y*q.z - z*q.y,
      y = + w*q.y - x*q.z + y*q.w + z*q.x,
      z = + w*q.z + x*q.y - y*q.x + z*q.w
    )
  
  def /(q: Quaternion): Quaternion = this * q.inverse
  
  def +(q: Quaternion): Quaternion = Quaternion(w+q.w, x+q.x, y+q.y, z+q.z)
  def -(q: Quaternion): Quaternion = Quaternion(w-q.w, x-q.x, y-q.y, z-q.z)
  
  def toGLRotationMatrix = {
    val m00 = 1.0 - 2.0*y*y - 2.0*z*z
    val m10 = 2.0*(x*y + w*z)
    val m20 = 2.0*(x*z - w*y)
    
    val m01 = 2.0*(x*y - w*z)
    val m11 = 1.0 - 2.0*x*x - 2.0*z*z
    val m21 = 2.0*(y*z + w*x)
    
    val m02 = 2.0*(x*z + w*y)
    val m12 = 2.0*(y*z - w*x)
    val m22 = 1.0 - 2.0*x*x - 2.0*y*y
    
    // this matrix is "visually transposed"	
    Array(
      m00, m10, m20, 0,
      m01, m11, m21, 0,
      m02, m12, m22, 0,
      0,   0,   0,   1
    )
  }
  
  def rotate(r: Vec3): Vec3 = {
    val m00 = 1.0 - 2.0*y*y - 2.0*z*z
    val m10 = 2.0*(x*y + w*z)
    val m20 = 2.0*(x*z - w*y)
    
    val m01 = 2.0*(x*y - w*z)
    val m11 = 1.0 - 2.0*x*x - 2.0*z*z
    val m21 = 2.0*(y*z + w*x)
    
    val m02 = 2.0*(x*z + w*y)
    val m12 = 2.0*(y*z - w*x)
    val m22 = 1.0 - 2.0*x*x - 2.0*y*y
    
    Vec3(m00*r.x + m01*r.y + m02*r.z,
         m10*r.x + m11*r.y + m12*r.z,
         m20*r.x + m21*r.y + m22*r.z)
  }
}
