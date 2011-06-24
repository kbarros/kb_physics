package kip.md

import kip.math.mutable


case class Tag(
  inter1: Seq[Interaction1] = Seq(),
  inter2: Seq[Interaction2] = Seq(),
  inter3: Seq[Interaction3] = Seq()
)


class Atom(var idx: Int, var tag: Tag, var mass: Double = 1d,
           val pos: mutable.Vec3 = mutable.Vec3.zero) {
  var wx, wy, wz: Int = 0 // wrapping indices
  val v = mutable.Vec3.zero
  var f = mutable.Vec3.zero // instantaneous force
  
  override def toString() = "Atom(%d, %s)".format(idx, pos)
}

