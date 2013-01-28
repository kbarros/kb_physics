package kip.projects.kondo_rkky

import java.awt.Color.ORANGE
import java.awt.Color.RED

import javax.swing.SwingUtilities
import kip.graphics.Bounds3d
import kip.graphics.RetainedScene
import kip.math.Vec3


class SpinViz(w: Int, h: Int) {
  val latDel1 = Vec3(0, 1, 0)
  val latDel3 = Vec3(1, 0, 0)

  val lat0 = Vec3(0, 0.5*math.sqrt(3)*(h-1), 0)
  val latN = lat0 + latDel1*(h-1) + latDel3*(w-1)
  
  val bds = Bounds3d(lat0, latN)

  def readSpin(x: Int, y: Int, field: Array[Double]): Vec3 = {
    require(x < w && y < w)
    val sx = field(0 + 3*(x + w*y))
    val sy = field(1 + 3*(x + w*y))
    val sz = field(2 + 3*(x + w*y))
    Vec3(sx, sy, sz)
  }
  
  def readPos(x: Int, y: Int): Vec3 = {
    lat0 + latDel3*x + latDel1*y
  }
  
  val spinPos: Array[Vec3] = {
    val ret = for (y <- (h-1) to 0 by -1;
         x <- 0 until w) yield {
      readPos(x, y)
    }
    ret.toArray
  }

  def spinDir(field: Array[Double]): Array[Vec3] = {
    val ret = for (y <- (h-1) to 0 by -1;
                   x <- 0 until w) yield {
      readSpin(x, y, field)
    }
    ret.toArray
  }

  def drawSpins(field: Array[Double], rs: RetainedScene) {

    val sd = spinDir(field)
    
    def omitSpin(i: Int) = {
      val x = i % w
      val y = i / w
      val frac = 0.25
      (x < frac*w || x+1 > (1-frac)*w || y < frac*w || y+1 > (1-frac)*w)
      false
    }
    
    val arrows = for (i <- 0 until w*w; if !omitSpin(i)) yield {
      val pos = spinPos(i) + Vec3(0, 0, 1)
      val spin = sd(i)
      val delta = spin*1.5
      val width = 0.3
      
      import java.awt.Color._
      val gray = new java.awt.Color(0, 0, 0, 50)

      new RetainedScene.Arrow(pos, delta, width, color1=ORANGE, color2=RED)
    }
    
    rs.bds = bds
    rs.drawables = Vector()
    rs.drawables ++= arrows
    
    SwingUtilities.invokeLater(new Runnable() {
      def run() = rs.display()
    })
  }
}
