package kip.md.apps

import java.awt._
import kip.md._
import kip.math.{Vec3, Quaternion}
import scala.math._

object LJTest {
  def main(args: Array[String]) {
    val L = 30
    val volume = new Volume.Cuboid(L, L, 0, periodic=true)
    val integrator = new Verlet(0.01)
    
    val pairlj = new LennardJones(eps=1, sigma_local=1, scaled_cutoff=3)
    val pairsoft = new PairSoft(eps=1, sigma_local=1)
    val tag = new Tag(inter2 = Seq(pairsoft))
    
    val rand = new util.Random(0)
    val atoms = Seq.tabulate(1000) { i:Int =>
      new Atom(idx=i, tag=tag, x=(L*rand.nextDouble), y=(L*rand.nextDouble), z=(L*rand.nextDouble))
    }
    
    // val atoms = Seq(
    //     new Atom(idx=0, tag=tag, x=0, y=1, z=0),
    //     new Atom(idx=1, tag=tag, x=0, y=(-1), z=0)
    // )
    
    val world = new World(volume, atoms, integrator)
    
    val viz = new Visualizer()
    
    println("energy: " + (world.potentialEnergy() + world.kineticEnergy()))
    
    for (i <- 0 until 100) {
      for (i <- 0 until 10)
        world.step()
      world.visualize(viz, radius = _ => 0.5, color = _ => java.awt.Color.BLUE)
    }
    
    kip.util.Util.time(for (i <- 0 until 500) world.step(), "500 steps")

    println("energy: " + (world.potentialEnergy() + world.kineticEnergy()))
  }
}
