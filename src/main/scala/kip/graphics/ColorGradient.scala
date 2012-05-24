package kip.graphics

import java.awt.Color
import math._

object ColorGradient {
  def yellowPurple(lo: Double, hi: Double) = {
    new ColorGradient(lo, hi,
        Array[(Double, (Int,Int,Int))] (
            (1-1.0,  (0, 0, 0)),
            (1-0.98, (10, 0, 50)),
            (1-0.95, (20, 0, 80)),
            (1-0.85, (61, 0, 130)), // blue
            (1-0.7,  (121, 20, 150)), // blue
            (1-0.5,  (190, 40, 90)), // solid red
            (1-0.35, (215, 90, 40)), // red
            (1-0.15, (235, 195, 80)), // yellow
            (1-0,    (255, 255, 255)))
    )
  }

  def blueRed(lo: Double, hi: Double) = {
    new ColorGradient(lo, hi,
        Array[(Double, (Int,Int,Int))] (
            (0,    (100, 100, 255)),
            (0.5,  (255, 255, 255)),
            (1,    (255, 100, 100)))
    )
  }
}


class ColorGradient(val lo: Double, val hi: Double, val controlPoints: Array[(Double, (Int,Int,Int))]) {
  val WHEEL_SIZE = 512
  val wheel = Array.tabulate(WHEEL_SIZE) { i =>
    val a = i.toDouble / WHEEL_SIZE
    val j = controlPoints.lastIndexWhere(_._1 <= a)
    val (v1, (r1, g1, b1)) = controlPoints(j)
    val (v2, (r2, g2, b2)) = controlPoints(j+1)
    
    val v = (a - v1) / (v2 - v1)
    val r = (r1*(1-v) + r2*v).toInt
    val g = (g1*(1-v) + g2*v).toInt
    val b = (b1*(1-v) + b2*v).toInt
    new Color(r, g, b)
  }
  
  
  def interpolate(x: Double) = {
    val v = (x - lo) / (hi - lo)
    val c = (WHEEL_SIZE*v).toInt
    wheel(min(max(c, 0), WHEEL_SIZE-1))
  }
}

