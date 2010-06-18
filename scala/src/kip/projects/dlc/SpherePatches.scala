package kip.projects.dlc

import scala.io.Source
import scala.annotation.tailrec
import java.io.File

import kip.math.Vec3
import kip.util.Util


/**
 * @param pos the center of patch
 * @param normal the surface normal
 * @param area the estimated voronoi area
 * @param curvature the mean curvature, 1/radius for the sphere
 */
case class Patch(pos: Vec3, normal: Vec3, area: Double, curvature: Double)


object SpherePatches {
  def fromFile(fname: String) = {
    println("Reading file "+fname)
    val lines = Source.fromFile(new File(fname)).getLines()
    val patches = lines.map(_.toDouble).grouped(3).map { case Seq(x,y,z) =>
      val pos = Vec3(x,y,z).normalize
      Patch(pos=pos, normal=pos, area=0, curvature=1)
    }
    new SpherePatches(patches.toArray, errorFraction=0.0)
  }
  
  def main(args: Array[String]) {
    val proj = "/home/kbarros/dev/projects/dielectric/packings"
    val n = 1472
    val typ = "maxvol"
    val p = fromFile(proj+"/"+n+"-"+typ+".txt").estimateArea(0.002)
    Util.writeStringToFile(p.toString, proj+"/geom_"+typ+"_"+n+".txt")
  }
}

class SpherePatches(val patches: Array[Patch], val errorFraction: Double) {
  val rand = new scala.util.Random(0)
  
  // given 'n' patches and 'm' samples, each patch is sampled approximately m/n times,
  // so approximate relative error is e = 1/sqrt(m/n) = sqrt(n/m)
  // this gives m = n / e^2
  def numSamplesForError(errorFraction: Double) = {
    if (errorFraction == 0) {
      0
    }
    else {
      val n = patches.size
      val e = errorFraction
      (n / (e*e)).toInt
    }
  }

  @tailrec
  final def randomPointOnSphere(): Vec3 = {
    def coord() = 2.0*rand.nextDouble() - 1.0
    val (x, y, z) = (coord(), coord(), coord())
    val d2 = x*x + y*y + z*z
    if (d2 < 1.0)
      Vec3(x, y, z).normalize
    else
      randomPointOnSphere()
  }
  
  def nearestPatchIndex(x: Vec3): Int = {
    var minDist2 = Double.MaxValue
    var minIdx = -1
    for (i <- patches.indices) {
      val p = patches(i)
      val dist2 = p.pos distance2 x
      if (dist2 < minDist2) {
        minDist2 = dist2
        minIdx = i
      }
    }
    minIdx
  }
  
  /**
   * Monte-Carlo estimate of each patch's Voronoi area
   * @param errorFraction target fractional error for each patch
   */
  def estimateArea(errorFraction: Double) = {
    val numSamples = numSamplesForError(errorFraction)
    val printFrequency = numSamples / 50
    val acc = Array.fill(patches.size)(0.0)
    
    Util.time ( {
      for (iter <- 0 until numSamples) {
        if (iter % printFrequency == 0)
          print(".")
        val x = randomPointOnSphere()
        val i = nearestPatchIndex(x)
        acc(i) += 1
      }
    }, "Evaluating %d samples".format(numSamples))
    
    val newPatches = Array.tabulate(patches.size) { i =>
      val sphereArea = 4.0*math.Pi
      val patchArea = (acc(i) / numSamples) * sphereArea
      patches(i).copy(area = patchArea)
    }
    new SpherePatches(newPatches, errorFraction)
  }
  
  override def toString = {
    val sb = new StringBuilder
    sb append
"""# Grid for sphere of radius 1 with %d samples and %g fractional error
# %d patches
# x y z nx ny nz area curvature
""".format(numSamplesForError(errorFraction), errorFraction, patches.size)
    for (p <- patches) {
      sb append "%14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f\n".format(
        p.pos.x, p.pos.y, p.pos.z,
        p.normal.x, p.normal.y, p.normal.z,
        p.area, p.curvature
      )
    }
    sb.toString
  }
}
