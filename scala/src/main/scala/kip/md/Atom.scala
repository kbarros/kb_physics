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

  // TODO: move this stuff into World, and accumulate force implicitly
  
  def potential1(world: World): Double = {
    var ret = 0.0
    for (i <- tag.inter1) {
      ret += i.potential(world, this)
    }
    ret
  }
  
  def potential2(world: World, that: Atom): Double = {
    var ret = 0.0
    for (i1 <- tag.inter2;
         i2 <- i1.compatibleInteractions(that.tag.inter2)) {
      ret += i1.potential(world, this, i2, that)
    }
    ret
  }
  
  
  def force1(world: World): Vec3 = {
    var fx, fy, fz = 0.0
    for (i1 <- tag.inter1) {
      val f = i1.force(world, this)
      fx += f.x
      fy += f.y
      fz += f.z
    }
    Vec3(fx, fy, fz)
  }

  def force2(world: World, that: Atom): (Vec3, Vec3) = {
    var f1x, f1y, f1z = 0.0
    var f2x, f2y, f2z = 0.0
    for (i1 <- tag.inter2;
         i2 <- i1.compatibleInteractions(that.tag.inter2)) {
      val (f1,f2) = i1.force(world, this, i2, that)
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

