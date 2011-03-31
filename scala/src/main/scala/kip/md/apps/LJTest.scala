package kip.md.apps

import java.awt._
import kip.md._
import kip.math.{Vec3, Quaternion}


object LJTest {
  def main(args: Array[String]) {
    val world = new World()
    
    val viz = new Visualizer()
    viz.setParticles(Seq(Visualizer.Sphere(Vec3(0.2, 0.8, 1.0), 0.02, Color.RED)))
//    Interpreter.start(("molviz", viz))
  }
}
