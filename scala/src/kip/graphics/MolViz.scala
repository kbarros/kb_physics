package kip.graphics

import kip.util.{Snapshot, LammpsParser}
import kip.util.Util.time
import kip.math.{Vec3}
import java.awt.{Color, Frame, FileDialog}
import java.awt.event.{WindowAdapter, WindowEvent}


object MolViz {
  def main(args: Array[String]) {
//    load()
    interpreter
  }
  
  def load(): MolViz = {
    val frame = new Frame("Select data file...")
    val d = new FileDialog(frame)
    println(d.setDirectory("/Users/kbarros/dev/repo/projects/dielectric/"))
    d.setVisible(true)
    var file = d.getFile()
    if (file != null) {
      file = d.getDirectory() + file
      val snaps = time(LammpsParser.readLammpsDump(file), "Reading dump-file")
      time(new MolViz(snaps), "Creating MolViz")
    }
    else {
      null
    }
  }

  def interpreter {
    val output = new java.io.PrintWriter(new java.io.OutputStreamWriter(Console.out))
    val repl = new scala.tools.nsc.InterpreterLoop(None, output)
    val settings = new scala.tools.nsc.Settings
    settings.usejavacp.value = true
    settings.classpath.value = "/Users/kbarros/dev/repo/scala/target/scala_2.8.1.RC3/classes"
    repl.main(settings)
  }

}

class MolViz(snaps: Seq[Snapshot]) {
  val colors = collection.mutable.Map[Int,Color]()
  val radii  = collection.mutable.Map[Int, Double]()
  val (frame, scene) = createGui()

  def goto(i: Int) {
    idx = i
    scene.triggerRepaint()
  }
  
  private var idx: Int = 0
  private def createGui(): (Frame, Scene) = {
    val scene = new Scene() {
      def drawContent(gfx: GfxGL) {
        println(java.lang.Thread.currentThread)
        val snap = snaps(idx)
        val bds = Bounds3d(snap.lo, snap.hi) // Box size depends on idx
        gfx.perspective3d(bds, rotation)
        gfx.setColor(Color.GREEN)
        gfx.drawCuboid(bds)
        for (i <- 0 until snap.x.size) {
          val pos = Vec3(snap.x(i), snap.y(i), snap.z(i))
          val typ = snap.typ(i).toInt
          val rad = radii.getOrElse(typ, 1.0)
          gfx.setColor(colors.getOrElse(typ, Color.RED))
          gfx.drawSphere(pos, rad)
        }
        gfx.ortho2d(Bounds3d(Vec3(0,0,0), Vec3(1,1,1)))
        gfx.setColor(Color.RED)
        gfx.rasterString("index " + idx, 0, 0)
      }
    }
    
    val frame = new Frame("JOGL HelloWorld2")
    frame.add(scene.component)
    frame.setSize(300, 300)
    frame.setVisible(true)
    frame.addWindowListener(new WindowAdapter() {
      override def windowClosing(e: WindowEvent) {
	System.exit(0)
      }
    })
    (frame, scene)
  }
}
