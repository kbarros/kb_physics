package kip.md.apps

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
      initDislocs = Some(d5.size)
    
    dislocHistory.append2(frameData.time, allDislocs.toDouble / initDislocs.get)
    bulkDislocHistory.append2(frameData.time, bulkDislocs.toDouble / initDislocs.get)
    dislocPlot.registerLines("All", dislocHistory, java.awt.Color.BLACK)
    dislocPlot.registerLines("Bulk", bulkDislocHistory, java.awt.Color.RED)
  }
}
