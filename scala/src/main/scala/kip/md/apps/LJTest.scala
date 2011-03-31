package kip.md.apps

import java.awt._
import kip.md._
import kip.math.{Vec3, Quaternion}
import scala.math._

object LJTest {
  def main(args: Array[String]) {
    
    val viz = new Visualizer()
    viz.setParticles(Seq(Visualizer.Sphere(Vec3(0.2, 0.8, 1.0), 0.02, Color.RED)))

    val volume = new Volume.Cuboid(1, 2, 0)
    val integrator = new Verlet(0.1)
    
    val lj = new LennardJones(eps=1, sigma=1, sigma_cutoff=pow(2, 1./6))
    val tag = new Tag(inter2 = Seq(lj))
    val atoms = Seq(
        new Atom(idx=0, tag=tag, x=0, y=0, z=0),
        new Atom(idx=1, tag=tag, x=1, y=1, z=0)
    )
    
    val world = new World(volume, atoms, integrator)
    
    world.step()
    
//    Interpreter.start(("molviz", viz))


    // val pair = new LennardJones()
    // val tag = new Tag(inter2 = Seq(pair))
    // val atoms = Seq(
    //   new Atom(idx=0, tag=tag, x=0, y=0, z=0),
    //   new Atom(idx=1, tag=tag, x=1, y=1, z=0)
    // )
    
    // val world = new World()
    // world.setSize(10, 10)
    // world.atoms = atoms
    
    // val viz = new Visualizer()
    
    // world.step()
    // world.visualize(viz, radius = _ => 1, color = _ => java.awt.Color.BLUE)
    
    // kip.util.Interpreter.start("world"->world, "viz"->viz)

  }
}
