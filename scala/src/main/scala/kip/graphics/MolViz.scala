package kip.graphics

import kip.util.{Snapshot, LammpsParser}
import kip.util.Util.time
import kip.util.Interpreter
import kip.math.{Vec3}

import java.awt.{BorderLayout, Color, Frame, FileDialog}
import java.awt.event.{WindowAdapter, WindowEvent}

import javax.swing.{JPanel, JSlider}
import javax.swing.event.{ChangeEvent, ChangeListener}


object RenderProperties {
  object ColorGradient {
    var colors = Array[(Double, Color)](
      0.0d->new Color(1f,0.2f,0f),
      0.1d->new Color(1f,0.4f,0f),
      0.5d->new Color(0.7f, 0.7f, 0.7f),
      0.9d->new Color(0f,0.4f,1f),
      1.0d->new Color(0f,0.2f,1f)
    )
    
    def apply(x: Double): Color = {
      def avgColors(c1: Color, c2: Color, x: Double): Color = {
        new Color((c1.getRed()*(1-x) + c2.getRed()*x).toFloat/255f,
                  (c1.getGreen()*(1-x) + c2.getGreen()*x).toFloat/255f,
                  (c1.getBlue()*(1-x) + c2.getBlue()*x).toFloat/255f)
      }
      
      for ((x0, c0) <- colors.headOption)
        if (x < x0) return c0
      
      for ((xn, cn) <- colors.lastOption)
        if (xn < x) return cn
      
      for (i <- 0 until colors.size-1) {
        val (x1, c1) = colors(i)
        val (x2, c2) = colors(i+1)
        if (x1 <= x && x <= x2) {
          return avgColors(c1, c2, (x - x1)/(x2-x1))
        }
      }
      
      Color.black
    }
  }

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
        ColorGradient(x)
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

  val spherepoint = new RenderProperties {
    def color(snap: Snapshot, atom: Int) = {
      val q = snap.q(atom)
      val qmin = -0.0012
      val qmax = 0.0012
        
      // map charge linearly to range [0, 1]
      val x = (q - qmin) / (qmax - qmin)
      
      // color gradient between red and blue
      ColorGradient(x)
    }
    
    def radius(snap: Snapshot, atom: Int) = {
      snap.typ(atom).toInt match {
        case 1 => 0.0 + 0.1 // surface patch
        case 2 => 1.0 - 0.1 // core particles
        case 4 => 0.1       // salt coion
      }
    }
    
    def sphereRes(snap: Snapshot, atom: Int) = {
      snap.typ(atom).toInt match {
        case 1 => 2 // surface patch
        case 2 => 3 // core particles
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
      // (home+"/Desktop/dlc-data/n100_v0.05_qr1_b400_p372_k0.1-1/dump3.dat", 1)
      // (home+"/Desktop/dlc-data/n100_v0.05_qr1_b400_p372_k10/dump2-0.gz", 1000)
      // (home+"/Desktop/dlc-data/n100_v0.05_qr1_b400_p372_k10-1/dump3-0.gz", 20)
      (home+"/Desktop/dlc-data/n100_v0.05_qr1_b400_p372_k10-2/dump3-0.gz", 20)
    }
    Interpreter.start(("molviz", makeMolviz(file, readEvery)))
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
}

class MolViz(val snaps: Seq[Snapshot], render: RenderProperties) {
  var idx: Int = 0
  
  val scene = new Scene() {
    def drawContent(gfx: GfxGL) {
      val snap = snaps(idx)
      val bds = Bounds3d(snap.lo, snap.hi) // Box size depends on idx
      // val bds = Bounds3d(Vec3(-1,-1,-1), Vec3(3,3,3))
      gfx.perspective3d(bds, rotation, translation*(bds.hi-bds.lo).norm)
      gfx.setColor(Color.GREEN)
      gfx.drawCuboid(bds)
      for (i <- 0 until snap.x.size) {
        val pos = Vec3(snap.x(i), snap.y(i), snap.z(i))
        gfx.setColor(render.color(snap, i))
        gfx.drawSphere(pos, render.radius(snap, i), render.sphereRes(snap, i))
//        println("drawing %f %f %f r=%f".format(pos.x, pos.y, pos.z, render.radius(snap,i)))
      }
      // val n = 100
      // for (i <- 0 until n) {
      //   val pos = Vec3(2*i/n, 2, 2)
      //   gfx.setColor(ColorGradient(i/100f))
      //   gfx.drawSphere(pos, 0.1, 1)
      // }
      
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

  def animate(molviz: MolViz) {
    for (i <- 500 until 2000) {
      molviz.goto(i)
      val comp = molviz.scene.canvas
      val im = kip.graphics.Utilities.captureJComponentImage(comp, comp.getWidth(), comp.getHeight())
      kip.graphics.Utilities.writeImage(im, "/Users/kbarros/Desktop/images/im"+i+".png")
    }
  }
}
