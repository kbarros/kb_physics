package kip.math.mutable


object Vec3 {
  def zero = Vec3(0,0,0)
}

case class Vec3(var x: Double, var y: Double, var z: Double) extends kip.math.Vec3 {
  def +=(v: kip.math.Vec3) {
    x += v.x
    y += v.y
    z += v.z
  }

  def -=(v: kip.math.Vec3) {
    x -= v.x
    y -= v.y
    z -= v.z
  }

  def reset() {
    x = 0
    y = 0
    z = 0
  }
}
