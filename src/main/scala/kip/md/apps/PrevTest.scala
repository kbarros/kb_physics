package kip.md.apps

import java.awt._
import kip.md._
import kip.math.{Vec3, mutable}
import scala.math._


object PrevTest {
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
      new Atom(idx=i, tag=tag1, pos=mutable.Vec3(L*rand.nextDouble, L*rand.nextDouble, L*rand.nextDouble))
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
      writer.write("%f %f %f\n".format(a.pos.x, a.pos.y, a.pos.z))
    }
    writer.close();
  }
  
  
  
  def test2() {
    
    val rand = new util.Random(0)
    
    val thermoDamp = Verlet.ThermoLangevin(temp=0, damp=10, rand)
    val integrator = new Verlet(dt=0.02, thermostat=thermoDamp)
    
    val sizeAsymmetryInit = 1.1
    val rows1 = 20
    val cols1 = 40
    val natoms1 = rows1*cols1
    val rows2 = rows1
    val cols2 = (cols1/sizeAsymmetryInit).toInt
    val natoms2 = rows2*cols2
    
    val r1 = pow(2, 1./6) / 2
    val r2 = r1 * sizeAsymmetryInit
    val off = 10
    
    // val wall1norm = Vec3(1,0,0)
    // val wall2norm = Vec3(-1,0,0)
    // var wall1posInit = Vec3(off-r1,0,0)
    // var wall2posInit = Vec3(off+2*r1*cols1,0,0)

    val wall1norm = Vec3(0,1,0)
    val wall2norm = Vec3(0,-1,0)
    var wall1posInit = Vec3(0, off-r1,0)
    var wall2posInit = Vec3(0, off+sqrt(3)*(r1*rows1 + r2*rows2), 0)
    
    val atoms = IndexedSeq.tabulate(natoms1+natoms2) { i:Int =>
      new Atom(idx=i, tag=null)
    }
    
    val volume = new Volume.Cuboid(2*r1*cols1+2*off, sqrt(3)*(rows1*r1+rows2*r2)+2*off, 0, periodic=true)
    val world = new World(volume, atoms, integrator)
    val viz = new Visualizer()
    viz.scene.translation = Vec3(0, 0, 0.5) // zoom in for improved movie frames
    viz.setBounds(volume.bounds)
    
    def setInteractions(sizeAsymmetry: Double, wall1pos: Vec3, wall2pos: Vec3) = {
      val pairlj1 = new LennardJones(eps=1, sigma_local=1, scaled_cutoff=3)
      val pairlj2 = new LennardJones(eps=1, sigma_local=sizeAsymmetry, scaled_cutoff=3)
      
      val wallA1 = new LJWall(pos=wall1pos, normal=wall1norm, sigma=0.5*1)
      val wallB1 = new LJWall(pos=wall2pos, normal=wall1norm, sigma=0.5*1)
      val wallA2 = new LJWall(pos=wall1pos, normal=wall2norm, sigma=0.5*sizeAsymmetry)
      val wallB2 = new LJWall(pos=wall2pos, normal=wall2norm, sigma=0.5*sizeAsymmetry)
      
      val tag1 = new Tag(inter1 = Seq(wallA1, wallB1), inter2 = Seq(pairlj1))
      val tag2 = new Tag(inter1 = Seq(wallA2, wallB2), inter2 = Seq(pairlj2))
      
      for (a <- atoms) {
        a.tag = if (a.idx < natoms1) tag1 else tag2
      }
    }
    
    def initGrid(x0: Double, y0: Double, r: Double, rows: Int, cols: Int, atoms: Seq[Atom]) = {
      assert(atoms.size == rows*cols)
      for (i <- 0 until rows;
           j <- 0 until cols) {
        val a = atoms(i*cols + j)
        a.pos.x = x0 + (i%2)*r + (2*r)*j
        a.pos.y = y0 + (sqrt(3)*r)*i
        a.pos.z = 0
      }
    }
    
    initGrid(x0=off, y0=off,                    r=r1, rows=rows1, cols=cols1, atoms.view(0, natoms1))
    initGrid(x0=off, y0=off+(sqrt(3)*r1)*rows1, r=r2, rows=rows2, cols=cols2, atoms.view(natoms1, natoms1+natoms2))
    
    for (i <- 0 until 800) {
      val asymSpeed = 0.0 // 0.002
      val asymDelta = asymSpeed*world.time
      val sizeAsymmetry = sizeAsymmetryInit + asymDelta
      
      val wallSpeed = 0.02
      val wallDelta = wallSpeed*world.time
      val wall1pos = wall1posInit + wall1norm*wallDelta
      val wall2pos = wall2posInit + wall2norm*wallDelta
      
      setInteractions(sizeAsymmetry, wall1pos, wall2pos)
      
      world.step(20)
      
      val neighbors = LJTest.neighborList(atoms, volume, cutoff=1.9)
      
      viz.setString("Atoms=%d, Temp=%f".format(atoms.size, world.temperature()))
      
      viz.setParticles(atoms.indices.map { i =>
        val a = atoms(i)
        val r = if (i < natoms1) r1 else r2
        val c = neighbors(i).size match {
          case 5 => Color.ORANGE
          case 6 => if (i < natoms1) Color.BLUE else Color.RED
          case 7 => Color.GREEN
          case _ => Color.PINK
        }
        Visualizer.Sphere(a.pos, radius=0.8*r, color=c)
      })
      
      viz.setWalls(Seq(
        Visualizer.Wall(wall1norm, wall1pos, Color.ORANGE),
        Visualizer.Wall(wall2norm, wall2pos, Color.ORANGE)
      ))
      
      viz.display()
      
      // javax.imageio.ImageIO.write(viz.scene.captureImage(), "PNG", new java.io.File("imgs/foo%d.png".format(i)))
    }
  }

}
