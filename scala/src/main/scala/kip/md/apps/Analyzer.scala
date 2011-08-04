package kip.md.apps
import kip.md.Visualizer
import kip.math.Vec3
import kip.graphics.Bounds3d
import java.awt.Color

object Analyzer extends App {
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
 
  val tempPlot = new scikit.graphics.dim2.Plot("Temperature")
  scikit.util.Utilities.frame(tempPlot.getComponent(), tempPlot.getTitle())
  val tempHistory = new scikit.dataset.DynamicArray()
 
  val dislocPlot = new scikit.graphics.dim2.Plot("Dislocations")
  scikit.util.Utilities.frame(dislocPlot.getComponent(), dislocPlot.getTitle())
  val dislocHistory = new scikit.dataset.DynamicArray()
  val bulkDislocHistory = new scikit.dataset.DynamicArray()
  var initDislocs = None: Option[Int]

  val dumpFiles = {
    val dfs = dir.listFiles().filter(_.getName.matches("""frame\d+\.gz"""))
    dfs.sortBy(df => """\d+""".r.findFirstIn(df.getName).get.toInt)
  }
  for (df <- dumpFiles) {
    val frameData = kip.util.Util.readObjectGz[FrameData](df.getPath)
    
    forceHistory.append2(frameData.time, frameData.force)
    forcePlot.registerLines("Data", forceHistory, java.awt.Color.BLACK)
    
    tempHistory.append2(frameData.time, frameData.temp)
    tempPlot.registerLines("Data", tempHistory, java.awt.Color.BLACK)
    
    val (d5, d7) = frameData.dislocations
    val allDislocs = d5.size
    val bulkDislocs = ((d5 zip d7).filter { case (i, j) => (i < frameData.natoms1 == j < frameData.natoms1) }).size
    
    if (initDislocs.isEmpty)
      initDislocs = Some(math.max(d5.size, 1))
    
    dislocHistory.append2(frameData.time, allDislocs.toDouble / initDislocs.get)
    bulkDislocHistory.append2(frameData.time, bulkDislocs.toDouble / initDislocs.get)
    dislocPlot.registerLines("All", dislocHistory, java.awt.Color.BLACK)
    dislocPlot.registerLines("Bulk", bulkDislocHistory, java.awt.Color.RED)
  }
  
  ///////////////////////////////////
  // Arrows
  val viz = new Visualizer(sizew=1000, sizeh=500)
  viz.scene.translation = Vec3(0, 0, 1.1) // zoom in for improved movie frames
  val bounds = {
    val sizeAsymmetry = 1.1
    val r1 = math.pow(2, 1./6) / 2 // radius of smaller atom
    val r2 = r1 * sizeAsymmetry
    val layers = 4 
    val cols1 = 200
    val rows1 = 20 
    val rows2 = rows1
    val cols2 = (cols1/sizeAsymmetry).toInt
    val layerWidth = math.sqrt(3)*(r1*rows1 + r2*rows2)
    val off = cols1*r1
    Bounds3d(Vec3.zero, Vec3(2*r1*cols1+2*off, layers*layerWidth+2*off, 0))
  }
  viz.setBounds(bounds)
  val targetNumArrows = 1000
  val fd0 = kip.util.Util.readObjectGz[FrameData](dumpFiles.head.getPath)
  val fd1 = kip.util.Util.readObjectGz[FrameData](dumpFiles.last.getPath)
  val rand = new util.Random(0)
  val arrowIndices = Seq.fill(targetNumArrows)(rand.nextInt(fd0.x.size)).distinct
  viz.arrowHeadSize = 1.7
  viz.arrows = for (i <- arrowIndices) yield {
    val p0 = Vec3(fd0.x(i), fd0.y(i), 0)
    val p1 = Vec3(fd1.x(i), fd1.y(i), 0)
    val scale = 0.5
    val del = (p0 - p1) * scale
    Visualizer.Arrow(from=p1+del, to=p1, normal=Vec3(0,0,1), color=(if (i < fd0.natoms1) Color.BLUE else Color.RED))
  }
  viz.display()

}
