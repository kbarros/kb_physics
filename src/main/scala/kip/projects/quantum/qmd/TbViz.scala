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
import scikit.graphics.dim2.Plot
import kip.projects.quantum.kpm.KPMUtil
import kip.projects.quantum.kpm.EnergyScale
import scikit.dataset.PointSet

object TbViz {
  def readSnap(file: File): TbMD.Snap = {
    kip.util.JacksonWrapper.deserialize[TbMD.Snap](file.slurp)
  }
  
  def readSnaps(dir: File, readEvery: Int): Array[TbMD.Snap] = {
    require(dir.isDirectory(), "Cannot load directory %s".format(dir))  
    val files = dir.listFiles()
    val indices = files.map { f =>
      f.getName().takeWhile(_ != '.').toInt
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
        gfx.drawSphere(pos, 1*Units.angstrom, 2)
      }
      
      gfx.ortho2dPixels()
      gfx.setColor(Color.RED)
      val pixOffset = 4
      gfx.rasterString(s"index $idx/${snaps.size-1} time ${snap.time} fs", pixOffset, pixOffset)

    }
  }
  
  val frame = {
    val frame = new Frame("Atoms")
    frame.addWindowListener(new WindowAdapter() {
      override def windowClosing(e: WindowEvent) {
        System.exit(0)
      }
    })

    val slider = new JSlider(0, snaps.size-1, 0)
    slider.addChangeListener(new ChangeListener() {
      def stateChanged(e: ChangeEvent) {
        goto(slider.getValue)
        scene.display()
      }
    })
    
    val panel = new JPanel(new BorderLayout())
    panel.add(scene.component, BorderLayout.CENTER);
    panel.add(slider, BorderLayout.SOUTH);
    
    frame.add(panel)
    frame.setSize(300, 300)
    frame.setVisible(true)
    frame
  }
  
  lazy val plot = {
    val plot = new Plot("Integrated density")
    scikit.util.Utilities.frame(plot.getComponent(), plot.getTitle())
    plot
  }
  
  def goto(i: Int) {
    idx = i
    scene.display()
    
    // integrated density 
    if (true) {
      val snap = snaps(idx)
      val es = new EnergyScale(snap.energyScale._1, snap.energyScale._2)
      val gamma = KPMUtil.momentTransform(snap.moments, 10*snap.moments.size)
      val (xp, irho) = KPMUtil.integratedDensityFunction(gamma, es)
      val data = new PointSet(xp, irho)
      plot.registerLines("Integrated density", data, Color.BLACK)
      val muIdx = xp.indexWhere(_ > snap.mu)
      val muData = new PointSet(Array(snap.mu), Array(irho(muIdx)))
      plot.registerBars("Mu", muData, Color.RED)
    }
  }
  
  goto(0)
  
  def animate() {
    for (i <- 500 until 2000) {
      goto(i)
      // val im = kip.graphics.Utilities.captureJComponentImage(comp, comp.getWidth(), comp.getHeight())
      // val im = scene.captureImage()
      // javax.imageio.ImageIO.write(im, "PNG", new java.io.File("foo6.png"))

    }
  }
}