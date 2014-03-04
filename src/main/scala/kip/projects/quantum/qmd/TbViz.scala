package kip.projects.quantum.qmd

import kip.util.{Snapshot, LammpsParser}
import kip.util.Util.time
import kip.util.Interpreter
import kip.math.{Vec3}
import kip.graphics._
import kip.enrich._
import java.awt.{BorderLayout, Color, Frame, FileDialog}
import java.awt.event.{WindowAdapter, WindowEvent}
import javax.swing.{JPanel, JSlider, JComponent}
import javax.swing.event.{ChangeEvent, ChangeListener}
import java.io.File

object TbViz {
  def readSnap(file: File): TbMD.Snap = {
    kip.util.JacksonWrapper.deserialize[TbMD.Snap](file.slurp)
  }
  
  def readSnaps(dir: File, readEvery: Int): Array[TbMD.Snap] = {
    require(dir.isDirectory(), "Cannot load directory %s".format(dir))  
    val files = dir.listFiles()
    val indices = files.map { f =>
      val name = f.getName()
      val idx = name.takeWhile(_ != '.').toInt
      println(s"Name = $name, idx=$idx")
      idx
    }
    (files zip indices).sortBy(_._2).filter(_._2 % readEvery == 0).map(x => readSnap(x._1))
  }
  
  def main(args: Array[String]) {
    if (args.size != 2) {
      println("Usage: TbViz <dirname> <readevery>")
    }
    val dir = args(0)
    val readEvery = args(1).toInt
    val snaps = readSnaps(new java.io.File(dir+"/dump"), readEvery)
    new TbViz(snaps)
  }
}


class TbViz(val snaps: Seq[TbMD.Snap]) {
  var idx: Int = 0
  
  val scene = new Scene() {
    def drawContent(gfx: GfxGL) {
      val snap = snaps(idx)
      val bds = {
        val lo = Vec3(snap.bdsLo(0), snap.bdsLo(1), snap.bdsLo(2))
        val hi = Vec3(snap.bdsHi(0), snap.bdsHi(1), snap.bdsHi(2))
        Bounds3d(lo, hi)
      }
      gfx.perspective3d(bds, rotation, translation*(bds.hi-bds.lo).norm)
      gfx.setColor(Color.GREEN)
      gfx.drawCuboid(bds)
      for (i <- 0 until snap.x.size/3) {
        val pos = Vec3(snap.x(3*i+0), snap.x(3*i+1), snap.x(3*i+2))
        gfx.setColor(Color.BLUE)
//        gfx.drawSphere(pos, 1*Units.angstrom, 2)
      }
      
      
      
      val origin = Vec3(0, 0, 0)
      val delta = Vec3(1, 1, 1)
      val width = 0.5
      val color1 = java.awt.Color.blue
      val color2 = java.awt.Color.red
      val a0 = delta.normalize 
      val a1 = {
        val alignedz = math.abs(a0.z) > (1 - 1e-6) 
        ((if (alignedz) Vec3(0, 1, 0) else Vec3(0, 0, 1)) cross a0).normalize    
      }
      val a2 = (a0 cross a1).normalize

      val l0 = delta.norm * 0.8 // length from base to waist
      val l1 = delta.norm - l0  // length from waist to end
      
      val x0 = origin
      val x1 = origin + a0 * l0
      val x2 = origin + a0 * (l0 + l1)
      
      val w0 = 0.5 * width * 0.3 // width at base
      val w1 = 0.5 * width       // width at waist
      
      gfx.ortho2dPixels()
      gfx.setColor(Color.RED)
      val pixOffset = 4
      gfx.rasterString(s"index $idx/${snaps.size-1} time ${snap.time} fs", pixOffset, pixOffset)

    }
  }
  
  val frame = buildFrame()
  
  def buildSlider(): JSlider = {
    val slider = new JSlider(0, snaps.size-1, 0)
    slider.addChangeListener(new ChangeListener() {
      def stateChanged(e: ChangeEvent) {
        goto(slider.getValue)
        scene.display()
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
    scene.display()
  }
  
  def animate(molviz: MolViz) {
    for (i <- 500 until 2000) {
      molviz.goto(i)

      // val im = kip.graphics.Utilities.captureJComponentImage(comp, comp.getWidth(), comp.getHeight())
      val im = scene.captureImage()
      javax.imageio.ImageIO.write(im, "PNG", new java.io.File("foo6.png"))

    }
  }
}