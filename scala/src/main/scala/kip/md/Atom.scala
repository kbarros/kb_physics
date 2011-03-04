package kip.md

import kip.math.Vec3


case class Tag(
  inter1: Seq[Interaction1] = Seq(),
  inter2: Seq[Interaction2] = Seq(),
  inter3: Seq[Interaction3] = Seq()
)


abstract class Atom extends PointGrid2d.Pt {
  val idx: Int
  var mass: Double
  var x, y, z: Double
  var wx, wy, wz: Int // wrapping indices
  var vx, vy, vz: Double // velocity
  var fx, fy, fz: Double // instantaneous force
  var tag: Tag
  
  def potential1: Double = {
    var ret = 0.0
    for (i <- tag.inter1) {
      ret += i.potential(this)
    }
    ret
  }
  
  def potential2(that: Atom): Double = {
    var ret = 0.0
    for (i1 <- tag.inter2;
         i2 <- i1.compatibleInteractions(that.tag.inter2)) {
      ret += i1.potential(this, i2, that)
    }
    ret
  }
  
  
  def force1: Vec3 = {
    var fx, fy, fz = 0.0
    for (i1 <- tag.inter1) {
      val f = i1.force(this)
      fx += f.x
      fy += f.y
      fz += f.z
    }
    Vec3(fx, fy, fz)
  }

  def force2(that: Atom): (Vec3, Vec3) = {
    var f1x, f1y, f1z = 0.0
    var f2x, f2y, f2z = 0.0
    for (i1 <- tag.inter2;
         i2 <- i1.compatibleInteractions(that.tag.inter2)) {
      val (f1,f2) = i1.force(this, i2, that)
      f1x += f1.x
      f1y += f1.y
      f1z += f1.z
      f2x += f2.x
      f2y += f2.y
      f2z += f2.z
    }
    (Vec3(f1x,f1y,f1z), Vec3(f2x,f2y,f2z))
  }
}

