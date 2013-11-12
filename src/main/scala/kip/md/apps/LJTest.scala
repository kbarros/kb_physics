package kip.md.apps

import java.awt._
import kip.md._
import kip.math.{Vec3, mutable}
import kip.util.JacksonWrapper._
import scala.collection.mutable.ArrayBuffer
import scala.math._


object LJTest {

/*
{
"layers": 2
"rows1": 80
"cols1": 80
"sizeAsymmetry": 1.0
}
*/
  case class Params(layers: Int, rows1: Int, cols1: Int, sizeAsymmetry: Double)

  def main(args: Array[String]) {
    require(args.size == 1, "Must pass configuration json file")
    val cfg = new java.io.File(args(0))
    
    val source = scala.io.Source.fromFile(args(0))
    
    val params = deserialize[Params](source.getLines().mkString)
    
    test3(layers=params.layers,
          rows1=params.rows1,
          cols1=params.cols1,
          sizeAsymmetry=params.sizeAsymmetry)
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
  
  
  def analyzeDislocations(t: Double, disloData: DislocationAnalysis, dislocationPairs: IndexedSeq[(Atom, Atom)]) {
    val points = dislocationPairs.map { case (a1, a2) => (a1.pos + a2.pos) / 2 }
    def dist(i: Vec3, j: Vec3): Double = (i - j).norm
    disloData.accum(t, points, points, dist)
    println("num dislocations: "+dislocationPairs.size)
  }

  def analyzeSurfaces(surfacePairs: IndexedSeq[(Atom, Atom)],
                      neighbors: IndexedSeq[Set[Int]], atoms: IndexedSeq[Atom]): Double = {
    def avg(s: Traversable[Double]) = s.sum / s.size
    def sqr(x: Double) = x*x

    val (atoms1, atoms2) = surfacePairs.unzip
    val surfaceAtoms = (atoms1 ++ atoms2).distinct

/*
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
    
    val variances = meanVar.unzip._2
    avg(variances)
*/
    surfaceAtoms.size
  }
  
  def test3(layers: Int, rows1: Int, cols1: Int, sizeAsymmetry: Double) {
    
    // Make target directory
    val dirname = "rows=%d_asym=%s_layers=%d".format(rows1, sizeAsymmetry, layers)
    val dir = new java.io.File(dirname)
    if (dir.exists) {
      if (!dir.isDirectory()) {
        println("File '%s' exists but is not directory".format(dir))
        sys.exit()
      }
      if (dir.list.size > 0) {
        println("WARNING: Directory: '%s' is not empty".format(dirname))
      }
    }
    if (!dir.exists) {
      dir.mkdir()
    }
    
    val rand = new util.Random(0)
    
    val thermoDamp = Verlet.ThermoLangevin(temp=0, damp=10, rand)
    val integrator = new Verlet(dt=0.02, thermostat=thermoDamp)
    

    val r1 = pow(2, 1.0/6) / 2 // radius of smaller atom
    val r2 = r1 * sizeAsymmetry

    val rows2 = rows1
    val cols2 = (cols1/sizeAsymmetry).toInt
    
    val layerWidth = sqrt(3)*(r1*rows1 + r2*rows2)
    
    val natoms1 = layers*rows1*cols1
    val natoms2 = layers*rows2*cols2
    val atomsPerLayer = cols1*rows1 + cols2*rows2
    
    val off = cols1*r1
    
    // val wall1norm = Vec3(1,0,0)
    // val wall2norm = Vec3(-1,0,0)
    // var wall1posInit = Vec3(off-r1,0,0)
    // var wall2posInit = Vec3(off+2*r1*cols1,0,0)

    val wall1norm = Vec3(0,1,0)
    val wall2norm = Vec3(0,-1,0)
    val wall1posInit = Vec3(0, off-r1,0)
    val wall2posInit = Vec3(0, off+layers*layerWidth, 0)
    
    val wallSpeed = 0.02
    val maxTime = 1000
    
    val stepsPerIter = 40
    val dumpFrequency = 2
    
    
    val atoms = IndexedSeq.tabulate(natoms1+natoms2) { i:Int =>
      new Atom(idx=i, tag=null)
    }
    
    val volume = new Volume.Cuboid(2*r1*cols1+2*off, layers*layerWidth+2*off, 0, periodic=true)
    val world = new World(volume, atoms, integrator)
    val viz = new Visualizer(sizew=1000, sizeh=500)
    viz.scene.translation = Vec3(0, 0, 1.0) // zoom in for improved movie frames
    viz.setBounds(volume.bounds)

//    val viz2 = new Visualizer(sizew=1000, sizeh=500)
//    viz2.scene.translation = Vec3(0, 0, 1.0) // zoom in for improved movie frames
//    viz2.setBounds(volume.bounds)
//    val targetNumArrows = 300
//    val arrowIndices = Seq.fill(targetNumArrows)(rand.nextInt(atoms.size)).distinct
    
    val forcePlot = new scikit.graphics.dim2.Plot("Force")
//    val roughPlot = new scikit.graphics.dim2.Plot("Roughening")
    val disloPlot = new scikit.graphics.dim2.Plot("Dislocations")
    scikit.util.Utilities.frame(forcePlot.getComponent(), forcePlot.getTitle())
//    scikit.util.Utilities.frame(roughPlot.getComponent(), roughPlot.getTitle())
    scikit.util.Utilities.frame(disloPlot.getComponent(), disloPlot.getTitle())
    val forceHistory = new scikit.dataset.DynamicArray()
//    val roughHistory = new scikit.dataset.DynamicArray()
    val disloHistory = new scikit.dataset.DynamicArray()
    val disloData = new DislocationAnalysis(tmin=0, tmax=maxTime, dt=100, rmin=0, rmax=50, dr=1, volume=1)
    
    def isType1(idx: Int): Boolean = idx < natoms1
    
    def setInteractions(sizeAsymmetry: Double, wall1pos: Vec3, wall2pos: Vec3): Seq[Seq[Interaction1]] = {
      val pairlj1 = new LennardJones(eps=1, sigma_local=1, scaled_cutoff=3)
      val pairlj2 = new LennardJones(eps=1, sigma_local=sizeAsymmetry, scaled_cutoff=3)
      
      val wallA1 = new LJWall(pos=wall1pos, normal=wall1norm, sigma=0.5*1)
      val wallB1 = new LJWall(pos=wall2pos, normal=wall2norm, sigma=0.5*1)
      val wallA2 = new LJWall(pos=wall1pos, normal=wall1norm, sigma=0.5*sizeAsymmetry)
      val wallB2 = new LJWall(pos=wall2pos, normal=wall2norm, sigma=0.5*sizeAsymmetry)
      
      val tag1 = new Tag(inter1 = Seq(wallA1, wallB1), inter2 = Seq(pairlj1))
      val tag2 = new Tag(inter1 = Seq(wallA2, wallB2), inter2 = Seq(pairlj2))
      
      for (a <- atoms) {
        a.tag = if (isType1(a.idx)) tag1 else tag2
      }
      
      Seq(Seq(wallA1, wallA2), Seq(wallB1, wallB2))
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
    val initialPos = atoms.map(_.pos.copy())
    
    def hasWrapped() = atoms.find(_.wx != 0).isDefined

    // =======================
    // Main loop
    //
    var iter = 0
    var netWallForce = Vec3.zero
    
    while (world.time < maxTime && !hasWrapped()) kip.util.Util.time("Iterating") {
      iter += 1
      
      // ----------------
      // World step
      //
      val wallDelta = wallSpeed*world.time
      val wall1pos = wall1posInit + wall1norm*wallDelta
      val wall2pos = wall2posInit + wall2norm*wallDelta
      val Seq(wallA, wallB) = setInteractions(sizeAsymmetry, wall1pos, wall2pos)
      world.step(stepsPerIter)
      val cutoff = 1.9
      val neighbors = neighborList(atoms, volume, cutoff)

      
      // ----------------
      // Force plot
      //
      val forceA = world.forceOnObject(wallA)
      val forceB = world.forceOnObject(wallB)
      val avgWallForce = (forceA.norm + forceB.norm) / 2.0
      netWallForce += forceA + forceB
//      println("force A = "+forceA)
//      println("force B = "+forceB)
//      println("avg force = "+avgWallForce)
//      println("time = "+world.time)
//      println("wall2pos = "+wall2pos)
      forceHistory.append2(world.time, avgWallForce)
      forcePlot.registerLines("Data", forceHistory, java.awt.Color.BLACK)
      
      
      // ----------------
      // Dislocations plot
      //
      
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
      // // Correlations between dislocations currently disabled
      // analyzeDislocations(world.time, disloData, dislocationPairs)
      // val disloResults = disloData.results
      // for (i <- disloResults.indices) {
      //   val (t, g) = disloResults(i)
      //   import java.awt.Color._
      //   val colors = Seq(BLACK, BLUE, RED, ORANGE, GRAY, YELLOW, GREEN)
      //   disloPlot.registerLines("Dislo "+t, new scikit.dataset.PointSet(g.elemCenters, g.elems), colors(i % colors.size))
      // }
      disloHistory.append2(world.time, dislocationPairs.size)
      disloPlot.registerLines("Dislo", disloHistory, java.awt.Color.BLACK)

      /*
      // ----------------
      // Roughening plot
      //
      val surfacePairs = {
        for (i <- atoms.indices;
             if (isType1(i));
             j <- neighbors(i);
             if (!isType1(j))) yield (atoms(i), atoms(j))
      }
      val roughening = analyzeSurfaces(surfacePairs, neighbors, atoms)
      roughHistory.append2(world.time, roughening)
      roughPlot.registerLines("Data", roughHistory, java.awt.Color.BLACK)
      */
      
      // ----------------
      // Atoms display
      //
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
      viz.display() // rendering can take 10 - 20% of the time
      javax.imageio.ImageIO.write(viz.scene.captureImage(), "PNG", new java.io.File(dirname+"/snap%d.png".format(iter)))
      
      
      // ----------------
      // Displacement display
      //
//      viz2.arrowHeadSize = 0.2
//      viz2.arrows = for (i <- arrowIndices) yield {
//        val p0 = initialPos(i)
//        val p1 = atoms(i).pos.copy()
//        val scale = 4
//        val del = (p0 - p1) / scale
//        Visualizer.Arrow(from=p1+del, to=p1, normal=Vec3(0,0,1), color=(if (isType1(i)) Color.BLUE else Color.RED))
//      }
//      viz2.display()

      
      // ----------------
      // Write frame data
      //
      if (iter % dumpFrequency == 0) {
        val filename = dirname+"/frame"+iter+".gz"
        val x = atoms.map(_.pos.x).toArray
        val y = atoms.map(_.pos.y).toArray
        val r = atoms.map(a => if (isType1(a.idx)) r1 else r2).toArray
        val (d5, d7) = dislocationPairs.unzip
        val frame = FrameData(index=iter, time=world.time, temp=world.temperature(), natoms1=natoms1, force=avgWallForce, x=x, y=y, r=r, dislocations=(d5.map(_.idx).toArray, d7.map(_.idx).toArray))
        kip.util.Util.writeObjectGz(filename, frame)
      }
    } // Simulation loop
    
  } // Simulation method
}
