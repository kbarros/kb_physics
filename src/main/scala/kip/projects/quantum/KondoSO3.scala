package kip.projects.quantum

import kip.graphics._
import kip.math.Vec3
import scala.collection.mutable.ArrayBuffer


class KondoSO3() {
  
  val tp = math.Pi
  val bds = Bounds3d(Vec3(-tp,-tp,-tp), Vec3(tp,tp,tp))
  val viz = new RetainedScene(bds, title="SO(3) path")
  
  // center and radius are in units of (2x2) spin plaquettes
  def loop(center: (Int, Int), radius: Int): Seq[(Int,Int)] = {
    // this (x,y) coordinates are in units of individual spin
    def mkPath(v0: (Int, Int), dv: (Int, Int), n: Int): Seq[(Int, Int)] = {
      val (x0, y0) = v0
//      println("Path corner at %d %d".format(x0, y0))
      val (dx, dy) = dv
      val ret = new ArrayBuffer[(Int, Int)]
      for (i <- 0 until n) {
        import KondoViz.{w,h}
        ret += (((x0+dx*i+w)%w, (y0+dy*i+h)%h))
      }
      ret
    }
    val (cx, cy) = center
    val l1 = mkPath((2*cx-2*radius, 2*cy-2*radius), (+2, 0), 2*radius)
    val l2 = mkPath((2*cx+2*radius, 2*cy-2*radius), (0, +2), 2*radius)
    val l3 = mkPath((2*cx+2*radius, 2*cy+2*radius), (-2, 0), 2*radius)
    val l4 = mkPath((2*cx-2*radius, 2*cy+2*radius), (0, -2), 2*radius)
    l1 ++ l2 ++ l3 ++ l4 :+ l1(0)
  }

  
  def matrixToAxisAngle(m: Array[Array[Double]]): Vec3 = {
    val epsilon = 0.01; // margin to allow for rounding errors
    val epsilon2 = 0.1; // margin to distinguish between 0 and 180 degrees
    
    if ((math.abs(m(0)(1)-m(1)(0))< epsilon)
      && (math.abs(m(0)(2)-m(2)(0))< epsilon)
      && (math.abs(m(1)(2)-m(2)(1))< epsilon)) {
        // singularity found
        // first check for identity matrix which must have +1 for all terms
        //  in leading diagonal and zero in other terms
        if ((math.abs(m(0)(1)+m(1)(0)) < epsilon2)
          && (math.abs(m(0)(2)+m(2)(0)) < epsilon2)
          && (math.abs(m(1)(2)+m(2)(1)) < epsilon2)
          && (math.abs(m(0)(0)+m(1)(1)+m(2)(2)-3) < epsilon2)) {
            // this singularity is identity matrix so angle = 0
            Vec3(0,0,0); // zero angle, arbitrary axis
        }
        else {
          // otherwise this singularity is angle = 180
          val angle = math.Pi;
          val xx = (m(0)(0)+1)/2;
          val yy = (m(1)(1)+1)/2;
          val zz = (m(2)(2)+1)/2;
          val xy = (m(0)(1)+m(1)(0))/4;
          val xz = (m(0)(2)+m(2)(0))/4;
          val yz = (m(1)(2)+m(2)(1))/4;
          if ((xx > yy) && (xx > zz)) { // m[0][0] is the largest diagonal term
            if (xx< epsilon) {
              Vec3(0, 0.7071, 0.7071)*angle
            }
            else {
              val x = math.sqrt(xx);
              val y = xy/x;
              val z = xz/x;
              Vec3(x, y, z)*angle
            }
          }
          else if (yy > zz) { // m[1][1] is the largest diagonal term
            if (yy< epsilon) {
              Vec3(0.7071, 0, 0.7071)*angle
            }
            else {
              val y = math.sqrt(yy);
              val x = xy/y;
              val z = yz/y;
              Vec3(x, y, z)*angle
            }   
          }
          else { // m[2][2] is the largest diagonal term so base result on this
            if (zz< epsilon) {
              Vec3(0.7071, 0.7071, 0)*angle
            }
            else {
              val z = math.sqrt(zz);
              val x = xz/z;
              val y = yz/z;
              Vec3(x, y, z)*angle
            }
          }
        }
    }
    else {
      // as we have reached here there are no singularities so we can handle normally
      var s = math.sqrt((m(2)(1) - m(1)(2))*(m(2)(1) - m(1)(2))
          +(m(0)(2) - m(2)(0))*(m(0)(2) - m(2)(0))
          +(m(1)(0) - m(0)(1))*(m(1)(0) - m(0)(1))); // used to normalise
      if (math.abs(s) < 0.001) s=1; 
      // prevent divide by zero, should not happen if matrix is orthogonal and should be
      // caught by singularity test above, but I've left it in just in case
      val angle = math.acos(( m(0)(0) + m(1)(1) + m(2)(2) - 1)/2);
      val x = (m(2)(1) - m(1)(2))/s;
      val y = (m(0)(2) - m(2)(0))/s;
      val z = (m(1)(0) - m(0)(1))/s;
      Vec3(x,y,z)*angle
    }
  }
  
  
  // rotation of unit cell
  //       (x+1, y+1)
  //   o - o
  //   | / |
  //   o - o
  // (x,y)
  // 
  // relative to S_{x,y}     = e_x
  //             S_{(x+1),y} = e_y
  //
  def siteToSO3(x: Int, y: Int, field: Array[R]): Vec3 = {
    val s1 = KondoViz.readSpin(x, y, field)
    val s2 = KondoViz.readSpin((x+1)%KondoViz.w, y, field)
    val s3 = (s1 cross s2).normalize
    val s4 = (s1 cross s3).normalize
    
    val mat = Array(s1.x, s1.y, s1.z,
                    s3.x, s3.y, s3.z,
                    s4.x, s4.y, s4.z)
    
    matrixToAxisAngle(mat.grouped(3).toArray)
  }
  
  
  def drawPath(xys: Array[(Int, Int)], field: Array[R]) {
    val pts = for ((x, y) <- xys) yield {
      siteToSO3(x, y, field)
    }
    
//    for (p <- pts) { println("%f %f %f".format(p.x, p.y, p.z)) }
    
    import java.awt.Color.{RED, GRAY}
    val colors = for (i <- pts.indices.dropRight(1)) yield {
      if (pts(i).distance(pts(i+1)) < math.Pi)
        RED
      else 
        GRAY
    }
    
    viz.drawables = Vector(new RetainedScene.Cuboid(bds), new RetainedScene.LineStrip(pts, colors.toArray))
    viz.display()
  }
}
