package kip.md.apps

import kip.md._
import kip.math.Vec3
import kip.graphics.Bounds3d
import java.awt.Color

object AnalyzerTemp extends App {
  if (args.size == 0) {
    println("Must pass a directory")
    sys.exit()
  }
  
  val dir = new java.io.File(args(0))
  if (!dir.exists || !dir.isDirectory()) {
    println("Path '%s' is not a directory".format(dir))
    sys.exit()
  }

  val forcePlot = new scikit.graphics.dim2.Plot("Force")
  scikit.util.Utilities.frame(forcePlot.getComponent(), forcePlot.getTitle())
  val forceHistory = new scikit.dataset.DynamicArray()

  val dumpFiles = {
    val dfs = dir.listFiles().filter(_.getName.matches("""frame\d+\.gz"""))
    dfs.sortBy(df => """\d+""".r.findFirstIn(df.getName).get.toInt)
  }
  
  val f0 = kip.util.Util.readObjectGz[FrameData](dumpFiles(0).getPath)
  
  // Set manually
  val layers: Int = 8
  val rows1: Int = 10
  val cols1: Int = 200
  val sizeAsymmetry: Double = 1.0
  
  val r1 = math.pow(2, 1./6) / 2 // radius of smaller atom
  val r2 = r1 * sizeAsymmetry
  val rows2 = rows1
  val cols2 = (cols1/sizeAsymmetry).toInt
  val layerWidth = math.sqrt(3)*(r1*rows1 + r2*rows2)
  val natoms1 = layers*rows1*cols1
  val natoms2 = layers*rows2*cols2
  require(natoms1 == f0.natoms1)
  require(natoms2 == f0.x.size - f0.natoms1)
  val atomsPerLayer = cols1*rows1 + cols2*rows2
  val off = cols1*r1  
  val wall1norm = Vec3(0,1,0)
  val wall2norm = Vec3(0,-1,0)
  val wall1posInit = Vec3(0, off-r1,0)
  val wall2posInit = Vec3(0, off+layers*layerWidth, 0)
  val wallSpeed = 0.02
  val atoms = IndexedSeq.tabulate(natoms1+natoms2) { i:Int =>
    new Atom(idx=i, tag=null)
  }
  val volume = new Volume.Cuboid(2*r1*cols1+2*off, layers*layerWidth+2*off, 0, periodic=true)
  val world = new World(volume, atoms, integrator=null)
  def isType1(idx: Int): Boolean = idx < natoms1
  
  def setInteractions(sizeAsymmetry: Double, wall1pos: Vec3, wall2pos: Vec3): Seq[Seq[Interaction1]] = {
    val wallA1 = new LJWall(pos=wall1pos, normal=wall1norm, sigma=0.5*1)
    val wallB1 = new LJWall(pos=wall2pos, normal=wall2norm, sigma=0.5*1)
    val wallA2 = new LJWall(pos=wall1pos, normal=wall1norm, sigma=0.5*sizeAsymmetry)
    val wallB2 = new LJWall(pos=wall2pos, normal=wall2norm, sigma=0.5*sizeAsymmetry)
    val tag1 = new Tag(inter1 = Seq(wallA1, wallB1))
    val tag2 = new Tag(inter1 = Seq(wallA2, wallB2))
    for (a <- atoms) {
      a.tag = if (isType1(a.idx)) tag1 else tag2
    }
    Seq(Seq(wallA1, wallA2), Seq(wallB1, wallB2))
  }
  
  for (df <- dumpFiles) {
    val frameData = kip.util.Util.readObjectGz[FrameData](df.getPath)
    
    val wallDelta = wallSpeed*(frameData.time - 0.8)
    val wall1pos = wall1posInit + wall1norm*wallDelta
    val wall2pos = wall2posInit + wall2norm*wallDelta
    val Seq(wallA, wallB) = setInteractions(sizeAsymmetry, wall1pos, wall2pos)
    
    for ((a, i) <- atoms.zipWithIndex) {
      a.pos.x = frameData.x(i)
      a.pos.y = frameData.y(i)
      a.pos.z = 0
    }
    
    val forceA = world.forceOnObject(wallA)
    val forceB = world.forceOnObject(wallB)
    val avgForce = (forceA.norm + forceB.norm) / 2
    
//    println("force A = "+forceA)
//    println("force B = "+forceB)
//    println("avg force = "+avgForce)
//    println("time = "+frameData.time)
//    println("wall2pos = "+wall2pos)
//    println()
    
    forceHistory.append2(frameData.time, avgForce)
    forcePlot.registerLines("Data", forceHistory, java.awt.Color.BLACK)
   }
}
