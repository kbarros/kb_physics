package kip.md

import kip.math.Vec3


trait Pt {
  var x, y, z: Double
  def toVec = Vec3(x, y, z)
}

case class Tag(interactions: Seq[Interaction2])


abstract case class Atom(x: Double, y: Double, z: Double) extends Pt {
  var vx, vy, vz: Double
  var wx, wy, wz: Int // wrapping indices
  var tag: Tag
  
  def potential(that: Atom): Double = {
    val d = Vec3(x-that.x, y-that.y, z-that.z)
    var ret = 0.0
    for (i1 <- tag.interactions;
         i2 <- i1.compatibleInteractions(that.tag.interactions)) {
      ret += i1.potential(i2, this, that)
    }
    ret
  }
}

