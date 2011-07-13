package kip.md.apps

import java.awt._
import kip.md._
import kip.math.{Vec3, mutable}
import scala.math._

object LJTest {
    
    
  def main(args: Array[String]) {
    test3(args.headOption)
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
  
  
  def neighborList(atoms: IndexedSeq[Atom], volume: Volume, cutoff: Double): IndexedSeq[Set[Int]] = {
    import kip.util.QHull._
    
    val tris = {
      val pos = new Array[Double](atoms.size*2)
      for (i <- 0 until atoms.size) {
        pos(2*i+0) = atoms(i).pos.x
        pos(2*i+1) = atoms(i).pos.y
      }
      delaunay2d(Vertices2d(pos))
    }
    
    val links = Array.fill(atoms.size+1)(Set[Int]())
    for (t <- 0 until tris.size/3) {
      def addEdge(i: Int, j: Int) {
        if (volume.distance2(atoms(i), atoms(j)) < cutoff*cutoff) {
          links(i) += j
          links(j) += i
        }
      }
      val a0 = tris(3*t+0)
      val a1 = tris(3*t+1)
      val a2 = tris(3*t+2)
      addEdge(a0,a1)
      addEdge(a1,a2)
      addEdge(a2,a0)
    }

    links
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
      
      val neighbors = neighborList(atoms, volume, cutoff=1.9)
      
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
  
  def analyzeDislocations(dislocationPairs: IndexedSeq[(Atom, Atom)]) {
    val points = dislocationPairs.map { case (a1, a2) => (a1.pos + a2.pos) / 2 }
    def dist(i: Int, j: Int) = (points(i) - points(j))/2d
    
    println("num dislocations: "+dislocationPairs.size)
  }

  def analyzeSurfaces(surfacePairs: IndexedSeq[(Atom, Atom)],
                      neighbors: IndexedSeq[Set[Int]], atoms: IndexedSeq[Atom]): Double = {
    def avg(s: Traversable[Double]) = s.sum / s.size
    def sqr(x: Double) = x*x

    val (atoms1, atoms2) = surfacePairs.unzip
    val surfaceAtoms = (atoms1 ++ atoms2).distinct
    val grouper = new Grouper(surfaceAtoms)
    for ((a1, a2) <- surfacePairs) grouper.bond(a1, a2)
    
    val meanVar: Seq[(Double, Double)] = for (group: Seq[Atom] <- grouper.groupMap.values.toSeq) yield {
      val type1 = group diff atoms2
      val meanY = avg(type1.map(_.pos.y))
      val varY  = avg(type1.map(a => sqr(meanY - a.pos.y)))
      (meanY, varY)
    }
//    val groups = grouper.groupMap.values.toSeq
//    println("num surfaces %d with num atoms %d".format(groups.size, groups(0).distinct.size))
    assert(meanVar.size == 1)
    
    val variances = meanVar.unzip._2
    avg(variances)
  }
  
  def test3(dirname: Option[String]) {
    
    val rand = new util.Random(0)
    
    val thermoDamp = Verlet.ThermoLangevin(temp=0, damp=10, rand)
    val integrator = new Verlet(dt=0.02, thermostat=thermoDamp)
    

    val sizeAsymmetryInit = 1.2
    val r1 = pow(2, 1./6) / 2
    val r2 = r1 * sizeAsymmetryInit

    val rows1 = 10
    val cols1 = 20 // 80
    val rows2 = rows1
    val cols2 = (cols1/sizeAsymmetryInit).toInt
    
    val layers = 1 // 4
    val layerWidth = sqrt(3)*(r1*rows1 + r2*rows2)
    
    val natoms1 = layers*rows1*cols1
    val natoms2 = layers*rows2*cols2
    val atomsPerLayer = cols1*rows1 + cols2*rows2
    
    val off = cols1*r1/2
    
    // val wall1norm = Vec3(1,0,0)
    // val wall2norm = Vec3(-1,0,0)
    // var wall1posInit = Vec3(off-r1,0,0)
    // var wall2posInit = Vec3(off+2*r1*cols1,0,0)

    val wall1norm = Vec3(0,1,0)
    val wall2norm = Vec3(0,-1,0)
    var wall1posInit = Vec3(0, off-r1,0)
    var wall2posInit = Vec3(0, off+layers*layerWidth, 0)
    
    val atoms = IndexedSeq.tabulate(natoms1+natoms2) { i:Int =>
      new Atom(idx=i, tag=null)
    }
    
    val volume = new Volume.Cuboid(2*r1*cols1+2*off, layers*layerWidth+2*off, 0, periodic=true)
    val world = new World(volume, atoms, integrator)
    val viz = new Visualizer()
    viz.scene.translation = Vec3(0, 0, 0.5) // zoom in for improved movie frames
    viz.setBounds(volume.bounds)

    val forcePlot = new scikit.graphics.dim2.Plot("Force")
    val roughPlot = new scikit.graphics.dim2.Plot("Roughening")
    scikit.util.Utilities.frame(forcePlot.getComponent(), forcePlot.getTitle())
    scikit.util.Utilities.frame(roughPlot.getComponent(), roughPlot.getTitle())
    val forceHistory = new scikit.dataset.DynamicArray()
    val roughHistory = new scikit.dataset.DynamicArray()
    
    def isType1(idx: Int): Boolean = idx < natoms1
    
    def setInteractions(sizeAsymmetry: Double, wall1pos: Vec3, wall2pos: Vec3) {
      val pairlj1 = new LennardJones(eps=1, sigma_local=1, scaled_cutoff=3)
      val pairlj2 = new LennardJones(eps=1, sigma_local=sizeAsymmetry, scaled_cutoff=3)
      
      val wallA1 = new LJWall(pos=wall1pos, normal=wall1norm, sigma=0.5*1)
      val wallB1 = new LJWall(pos=wall2pos, normal=wall1norm, sigma=0.5*1)
      val wallA2 = new LJWall(pos=wall1pos, normal=wall2norm, sigma=0.5*sizeAsymmetry)
      val wallB2 = new LJWall(pos=wall2pos, normal=wall2norm, sigma=0.5*sizeAsymmetry)
      
      val tag1 = new Tag(inter1 = Seq(wallA1, wallB1), inter2 = Seq(pairlj1))
      val tag2 = new Tag(inter1 = Seq(wallA2, wallB2), inter2 = Seq(pairlj2))
      
      for (a <- atoms) {
        a.tag = if (isType1(a.idx)) tag1 else tag2
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
    
    for (layer <- 0 until layers) {
      initGrid(x0=off, y0=off+layer*layerWidth,
               r=r1, rows=rows1, cols=cols1, atoms.view(layer*rows1*cols1, (layer+1)*rows1*cols1))
    }
    for (layer <- 0 until layers) {
      initGrid(x0=off, y0=off+layer*layerWidth+(sqrt(3)*r1)*rows1,
               r=r2, rows=rows2, cols=cols2, atoms.view(natoms1+layer*rows2*cols2, natoms1+(layer+1)*rows2*cols2))
    }
    
    def averageWallForce(): Double = {
      val inter1s = atoms.flatMap(_.tag.inter1).distinct
      val force1s = inter1s.map(inter => (world.forceOnObject(inter)).norm)
      val wallForce = force1s.sum / force1s.size
      wallForce
    }
    
    for (i <- 0 until 1200) kip.util.Util.time2("Iterating") {
      val wallSpeed = 0.02
      val wallDelta = wallSpeed*world.time
      val wall1pos = wall1posInit + wall1norm*wallDelta
      val wall2pos = wall2posInit + wall2norm*wallDelta
      setInteractions(sizeAsymmetryInit, wall1pos, wall2pos)
      
      world.step(20)

      forceHistory.append2(world.time, averageWallForce())
      forcePlot.registerLines("Data", forceHistory, java.awt.Color.BLACK)
      
      val cutoff = 1.9
      val neighbors = neighborList(atoms, volume, cutoff)

      val dislocationPairs = {
        for (i <- atoms.indices;
             if (neighbors(i).size == 7);
             neigh5 = neighbors(i).filter(neighbors(_).size == 5)
             if (neigh5.nonEmpty)) yield {
          val dists = neigh5.map(j => volume.distance2(atoms(i), atoms(j)))
          val (jmin, dist) = (neigh5 zip dists).minBy(_._2)
          (atoms(i), atoms(jmin))
        }
      }
      analyzeDislocations(dislocationPairs)
      
      val surfacePairs = {
        for (i <- atoms.indices;
             if (isType1(i));
             j <- neighbors(i);
             if (!isType1(j))) yield (atoms(i), atoms(j))
      }
      val roughening = analyzeSurfaces(surfacePairs, neighbors, atoms)
      roughHistory.append2(world.time, roughening)
      roughPlot.registerLines("Data", roughHistory, java.awt.Color.BLACK)
      
      viz.setString("Time=%d, Atoms=%d, Temp=%f".format(world.time.toInt, atoms.size, world.temperature()))
      
      viz.setParticles(atoms.indices.map { i =>
        val a = atoms(i)
        val r = if (isType1(i)) r1 else r2
        val c = neighbors(i).size match {
          case 5 => Color.ORANGE
          case 6 => if (isType1(i)) Color.BLUE else Color.RED
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

      // dirname.foreach { dirname => 
      //   javax.imageio.ImageIO.write(viz.scene.captureImage(), "PNG", new java.io.File(dirname+"/snap%d.png".format(i)))
      // }
    }
  }

}
