package kip.graphics

import kip.util.{Snapshot, LammpsParser}
import kip.util.Util.time
import kip.math.{Vec3}
import java.awt.{Color, Frame, FileDialog}
import java.awt.event.{WindowAdapter, WindowEvent}



object Render {
  val basic = new Render {
    def color(snap: Snapshot, atom: Int) = Color.red
    def radius(snap: Snapshot, atom: Int) = 1.0
  }
  
  val nano = new Render {
    def color(snap: Snapshot, atom: Int) = {
      val q = snap.q(atom) 
      val qmin = -1
      val qmax = 1
      
      // interpolating value in [0,1]
      val x = ((q - qmin) / (qmax - qmin)).toFloat
      
      // color gradient between red and blue
      new Color(x, 0f, 1f-x)
    }
    
    def radius(snap: Snapshot, atom: Int) = {
      snap.typ(atom).toInt match {
        case 1 => 0.3 // surface patch
        case 2 => 3 // core particles
        case 3 => 4 // salt counterions
        case 4 => 0.5 // salt coion
      }
    }
  }
}

abstract class Render {
  def color(snap: Snapshot, atom: Int): Color
  def radius(snap: Snapshot, atom: Int): Double
}


object MolViz {
  def main(args: Array[String]) {
    // interpreter()
    val file = "/Users/kbarros/dev/repo/projects/dielectric/cylinder.L20/dump.dat"
    interpreter(("molviz", makeMolviz(file)))
  }
  
  def makeMolviz(file: String): MolViz = {
    val snaps = time(LammpsParser.readLammpsDump(file), "Reading '%s'".format(file))
    time(new MolViz(snaps, Render.nano), "Creating MolViz")
  }
  
  def load(): MolViz = {
    val frame = new Frame("Select data file...")
    val d = new FileDialog(frame)
    println(d.setDirectory("/Users/kbarros/dev/repo/projects/dielectric/"))
    d.setVisible(true)
    var file = d.getFile()
    if (file != null)
      makeMolviz(d.getDirectory() + file)
    else
      null
  }

  def interpreter(bindings: (String, Any)*) {
    val output = new java.io.PrintWriter(new java.io.OutputStreamWriter(Console.out))
    val repl = new scala.tools.nsc.InterpreterLoop(None, output)
    val settings = new scala.tools.nsc.Settings
    settings.usejavacp.value = true
    settings.classpath.value = "/Users/kbarros/dev/repo/scala/target/scala_2.8.1.RC3/classes"
    repl.settings = settings
    repl.createInterpreter()
    val varStr = bindings.unzip._1.mkString("[",",","]")
    time(bindings.foreach{ case (k,v) => repl.injectOne(k, v) }, "Binding values "+varStr)
    repl.in = scala.tools.nsc.interpreter.InteractiveReader.createDefault(repl.interpreter)
    try {
      // it is broken on startup; go ahead and exit
      if (repl.interpreter.reporter.hasErrors) return
      repl.printWelcome()
      repl.repl()
    } finally {
      repl.closeInterpreter()
    }
  }
}

class MolViz(snaps: Seq[Snapshot], render: Render) {
  var idx: Int = 0
  
  val (frame, scene) = {
    val scene = new Scene() {
      def drawContent(gfx: GfxGL) {
        val snap = snaps(idx)
        val bds = Bounds3d(snap.lo, snap.hi) // Box size depends on idx
        gfx.perspective3d(bds, rotation, translation*(bds.hi-bds.lo).norm)
        gfx.setColor(Color.GREEN)
        gfx.drawCuboid(bds)
        for (i <- 0 until snap.x.size) {
          val pos = Vec3(snap.x(i), snap.y(i), snap.z(i))
          gfx.setColor(render.color(snap, i))
          gfx.drawSphere(pos, render.radius(snap, i))
        }
        gfx.ortho2d(Bounds3d(Vec3(0,0,0), Vec3(1,1,1)))
        gfx.setColor(Color.RED)
        gfx.rasterString("index " + idx, 0, 0)
      }
    }
    
    val frame = new Frame("JOGL HelloWorld2")
    frame.addWindowListener(new WindowAdapter() {
      override def windowClosing(e: WindowEvent) {
	System.exit(0)
      }
    })
    frame.add(scene.component)
    frame.setSize(300, 300)
    frame.setVisible(true)
    (frame, scene)
  }
  
  def goto(i: Int) {
    idx = i
    scene.triggerRepaint()
  }
  
}
