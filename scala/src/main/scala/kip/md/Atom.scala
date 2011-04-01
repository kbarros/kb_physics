package kip.md

import kip.math.Vec3


case class Tag(
  inter1: Seq[Interaction1] = Seq(),
  inter2: Seq[Interaction2] = Seq(),
  inter3: Seq[Interaction3] = Seq()
)


class Atom(var idx: Int, var tag: Tag, var mass: Double = 1d,
           var x: Double = 0, var y: Double = 0, var z: Double = 0) extends PointGrid2d.Pt {
  var wx, wy, wz: Int = 0 // wrapping indices
  var vx, vy, vz: Double = 0 // velocity
  var fx, fy, fz: Double = 0 // instantaneous force
  
  override def toString() = "Atom(%d, %s)".format(idx, pos)
  
  def pos = Vec3(x, y, z)

}

