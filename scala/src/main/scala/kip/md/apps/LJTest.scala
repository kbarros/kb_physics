package kip.md.apps

import java.awt._
import kip.md._
import kip.math.{Vec3, Quaternion}
import scala.math._

object LJTest {
    
    
  def main(args: Array[String]) {
    test2()
  }

  def test1() { 
    val L = 30
    val volume = new Volume.Cuboid(L, L, 0, periodic=true)
    
    val rand = new util.Random(0)
    
    val thermoDamp = Verlet.ThermoLangevin(temp=0, damp=1, rand)
    val integrator = new Verlet(dt=0.01, thermostat=thermoDamp)
    
    val pairlj = new LennardJones(eps=1, sigma_local=1, scaled_cutoff=3)
    val pairsoft = new PairSoft(eps=1, sigma_local=1)
    val tag1 = new Tag(inter2 = Seq(pairsoft))
    val tag2 = new Tag(inter2 = Seq(pairlj))
    
    val atoms = Seq.tabulate(800) { i:Int =>
      new Atom(idx=i, tag=tag1, x=(L*rand.nextDouble), y=(L*rand.nextDouble), z=(L*rand.nextDouble))
    }
    
    val world = new World(volume, atoms, integrator)
    
    val viz = new Visualizer()
    viz.setBounds(volume.bounds)
    
    for (i <- 0 until 100) {
      world.step(10)
      viz.setString("Atoms=%d, Temp=%f".format(atoms.size, world.temperature()))
      viz.setParticles(atoms.map(a => Visualizer.Sphere(a.pos, radius=0.5, color=Color.BLUE)))
    }
    
    atoms.foreach {_.tag = tag2}
    integrator.thermostat=Verlet.ThermoLangevin(temp=1, damp=1, rand)
    
    for (i <- 0 until 1000) {
      world.step(10)
      viz.setString("Atoms=%d, Temp=%f".format(atoms.size, world.temperature()))
      viz.setParticles(atoms.map(a => Visualizer.Sphere(a.pos, radius=0.5, color=Color.BLUE)))
    }

    // kip.util.Util.time(for (i <- 0 until 500) world.step(), "500 steps")
  }
  
  
  def writeConfig(file: String, atoms: Seq[Atom]) {
    import java.io._
    val writer = new BufferedWriter(new FileWriter(file))
    writer.write(atoms.size+"\n")
    for (a <- atoms) {
      writer.write("%f %f %f\n".format(a.x, a.y, a.z))
    }
    writer.close();
  }
  
  
  def test2() {
    
    val rand = new util.Random(0)
    
    val thermoDamp = Verlet.ThermoLangevin(temp=0, damp=10, rand)
    val integrator = new Verlet(dt=0.02, thermostat=thermoDamp)
    
    val sizeAsymmetry = 1.0
    val rows1 = 20
    val cols1 = 40
    val natoms1 = rows1*cols1
    val rows2 = rows1
    val cols2 = (cols1/sizeAsymmetry).toInt
    val natoms2 = rows2*cols2

    val atoms = IndexedSeq.tabulate(natoms1+natoms2) { i:Int =>
      new Atom(idx=i, tag=null)
    }
    
    def setInteractions(sizeAsymmetry: Double) = {
      val pairlj1 = new LennardJones(eps=1, sigma_local=1, scaled_cutoff=3)
      val pairlj2 = new LennardJones(eps=1, sigma_local=sizeAsymmetry, scaled_cutoff=3)
      val tag1 = new Tag(inter2 = Seq(pairlj1))
      val tag2 = new Tag(inter2 = Seq(pairlj2))
      for (a <- atoms) {
        a.tag = if (a.idx < natoms1) tag1 else tag2
      }
    }
    
    def initGrid(x0: Double, y0: Double, r: Double, rows: Int, cols: Int, atoms: Seq[Atom]) = {
      assert(atoms.size == rows*cols)
      for (i <- 0 until rows;
           j <- 0 until cols) {
        val a = atoms(i*cols + j)
        a.x = x0 + (i%2)*r + (2*r)*j
        a.y = y0 + (sqrt(3)*r)*i
        a.z = 0
      }
    }
    
    setInteractions(sizeAsymmetry)
    
    val r1 = pow(2, 1./6) / 2
    val r2 = r1 * sizeAsymmetry
    val off = 10
    initGrid(x0=off, y0=off,                    r=r1, rows=rows1, cols=cols1, atoms.view(0, natoms1))
    initGrid(x0=off, y0=off+(sqrt(3)*r1)*rows1, r=r2, rows=rows2, cols=cols2, atoms.view(natoms1, natoms1+natoms2))
    
    val volume = new Volume.Cuboid(2*r1*cols1+2*off, sqrt(3)*(rows1*r1+rows2*r2)+2*off, 0, periodic=true)
    
    val world = new World(volume, atoms, integrator)
    
    val viz = new Visualizer()
    viz.scene.translation = Vec3(0, 0, 0.5)
    viz.setBounds(volume.bounds)
    
    for (i <- 0 until 500) {
      world.step(10)
      viz.setString("Atoms=%d, Temp=%f".format(atoms.size, world.temperature()))
      viz.setParticles(atoms.map { a =>
        val (c, r) = if (a.idx < natoms1) (Color.BLUE, r1) else (Color.RED, r2)
        Visualizer.Sphere(a.pos, radius=0.8*r, color=c)
      })
      viz.display()
      
      setInteractions((1+0.002*world.time)*sizeAsymmetry)
      
      javax.imageio.ImageIO.write(viz.scene.captureImage(), "PNG", new java.io.File("imgs/foo%d.png".format(i)))
    }
  }
}
