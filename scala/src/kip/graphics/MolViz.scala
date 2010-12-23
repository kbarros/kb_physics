package kip.graphics

import kip.util.{Snapshot, LammpsParser}
import kip.util.Util.time
import kip.math.{Vec3}

import java.awt.{BorderLayout, Color, Frame, FileDialog}
import java.awt.event.{WindowAdapter, WindowEvent}

import javax.swing.{JPanel, JSlider}
import javax.swing.event.{ChangeEvent, ChangeListener}


object RenderProperties {
  val basic = new RenderProperties {
    def color(snap: Snapshot, atom: Int) = Color.red
    def radius(snap: Snapshot, atom: Int) = 1.0
    def sphereRes(snap: Snapshot, atom: Int) = 1
  }
  
  val nano = new RenderProperties {
    def color(snap: Snapshot, atom: Int) = {
      if (snap.q == null) {
        snap.typ(atom).toInt match {
          case 1 => null // surface patch
          case 2 => Color.blue // core particles
          case 3 => Color.red // salt counterions
          case 4 => Color.green // salt coion
        }
      }
      else {
        val q = snap.q(atom) 
        val qmin = -0.01
        val qmax = 0.01
        
        // map charge linearly to range [0, 1]
        val x = (q - qmin) / (qmax - qmin)
        
        // color gradient between red and blue
        def saturate(x: Double) = math.min(math.max(x, 0d), 1d).toFloat
        new Color(saturate(x), 0f, saturate(1-x))
      }
    }
    
    def radius(snap: Snapshot, atom: Int) = {
      snap.typ(atom).toInt match {
        case 1 => 0.0 + 0.5 // surface patch
        case 2 => 3.5 - 0.5 // core particles
        case 3 => 0.5 // salt counterions
        case 4 => 0.5 // salt coion
      }
    }
    
    def sphereRes(snap: Snapshot, atom: Int) = {
      snap.typ(atom).toInt match {
        case 1 => 2 // surface patch
        case 2 => 3 // core particles
        case 3 => 3 // salt counterions
        case 4 => 3 // salt coion
      }
    }
  }
}

abstract class RenderProperties {
  def color(snap: Snapshot, atom: Int): Color
  def radius(snap: Snapshot, atom: Int): Double
  def sphereRes(snap: Snapshot, atom: Int): Int
}


object MolViz {
  def main(args: Array[String]) {
    val (file, readEvery) = if (args.size > 0) {
      (args(0).toString, args(1).toInt)
    }
    else {
      val home = System.getProperty("user.home")
      // interpreter()
      // (home+"/dev/repo/projects/dielectric/cylinder.L20/dump.dat", 1)
      // (home+"/Desktop/dlc-data/n100_v0.05_qr1_b400_p372_k0.1/dump2-0.gz", 100)
      // (home+"/Desktop/dlc-data/n100_v0.05_qr1_b400_p0_k1/dump2-0.gz", 10)
      // (home+"/Desktop/dlc-data/n100_v0.05_qr1_b400_p372_k10/dump2-0.gz", 100)
      (home+"/Desktop/dlc-data/n100_v0.05_qr1_b400_p372_k0.1-1/dump3.dat", 1)
    }
    interpreter(("molviz", makeMolviz(file, readEvery)))
  }
  
  def makeMolviz(file: String, readEvery: Int): MolViz = {
    val snaps = time(LammpsParser.readLammpsDump(file, readEvery=readEvery), "Reading '%s'".format(file))
    time(new MolViz(snaps, RenderProperties.nano), "Creating MolViz")
  }
  
  def load(): MolViz = {
    val frame = new Frame("Select data file...")
    val d = new FileDialog(frame)
    println(d.setDirectory("/Users/kbarros/dev/repo/projects/dielectric/"))
    d.setVisible(true)
    var file = d.getFile()
    if (file != null)
      makeMolviz(d.getDirectory() + file, readEvery=1)
    else
      null
  }

  def interpreter(bindings: (String, Any)*) {
    val output = new java.io.PrintWriter(new java.io.OutputStreamWriter(Console.out))
    val repl = new scala.tools.nsc.InterpreterLoop(None, output)
    val settings = new scala.tools.nsc.Settings
    settings.usejavacp.value = true
    settings.classpath.value = System.getProperty("scala.class.path")
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

class MolViz(val snaps: Seq[Snapshot], render: RenderProperties) {
  var idx: Int = 0
  
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
        gfx.drawSphere(pos, render.radius(snap, i), render.sphereRes(snap, i))
      }
      gfx.ortho2dPixels()
      gfx.setColor(Color.RED)
      val pixOffset = 4
      gfx.rasterString("index %d/%d, time %d".format(idx, snaps.size-1, snap.time), pixOffset, pixOffset)
    }
  }
  
  val frame = buildFrame()
  
  def buildSlider(): JSlider = {
    val slider = new JSlider(0, snaps.size-1, 0)
    slider.addChangeListener(new ChangeListener() {
      def stateChanged(e: ChangeEvent) {
	goto(slider.getValue)
	scene.triggerRepaint()
      }
    })
    slider
  }
  
  def buildFrame(): Frame = {
    val frame = new Frame("Molecular Visualizer")
    frame.addWindowListener(new WindowAdapter() {
      override def windowClosing(e: WindowEvent) {
        System.exit(0)
      }
    })

    val panel = new JPanel(new BorderLayout())
    panel.add(scene.component, BorderLayout.CENTER);
    panel.add(buildSlider(), BorderLayout.SOUTH);
    
    frame.add(panel)
    frame.setSize(300, 300)
    frame.setVisible(true)
    frame
  }
  
  def goto(i: Int) {
    idx = i
    scene.triggerRepaint()
  }
}
