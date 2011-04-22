package kip.md

import java.awt.{Color, Frame, BorderLayout}
import javax.swing.JPanel
import java.awt.event._
import kip.graphics.{Scene, Bounds3d, GfxGL}
import kip.math.{Vec3, Quaternion}
import kip.util.Interpreter

object Visualizer {
  case class Sphere(pt: Vec3, radius: Double, color: Color, resolution: Int = 2)
}  

class Visualizer {
  
  var particles = Seq[Visualizer.Sphere]()
  var bds = Bounds3d(Vec3(0,0,0), Vec3(1,1,1))
  var rasterString = ""

  def setBounds(bds: Bounds3d) {
    this.bds = bds
  }
  
  def setParticles(particles: Seq[Visualizer.Sphere]) {
    this.particles = particles
  }
  
  def setString(str: String) {
    rasterString = str
  }
  
  def display() {
    scene.display()
  }
  
  val scene = new Scene() {
    def drawContent(gfx: GfxGL) {
      gfx.perspective3d(bds, rotation, translation*(bds.hi-bds.lo).norm)
      gfx.setColor(Color.GREEN)
      gfx.drawCuboid(bds)
      for (p <- particles) {
        gfx.setColor(p.color)
        gfx.drawSphere(p.pt, p.radius, p.resolution)
      }
      gfx.ortho2dPixels()
      gfx.setColor(Color.RED)
      val pixOffset = 4
      gfx.rasterString(rasterString, pixOffset, pixOffset)
    }
  }
  
  val frame = buildFrame()
  
  // def buildSlider(): JSlider = {
  //   val slider = new JSlider(0, snaps.size-1, 0)
  //   slider.addChangeListener(new ChangeListener() {
  //     def stateChanged(e: ChangeEvent) {
  //       goto(slider.getValue)
  //       scene.triggerRepaint()
  //     }
  //   })
  //   slider
  // }
  
  def buildFrame(): Frame = {
    val frame = new Frame("Molecular Visualizer")
    frame.addWindowListener(new WindowAdapter() {
      override def windowClosing(e: WindowEvent) {
        System.exit(0)
      }
    })

    val panel = new JPanel(new BorderLayout())
    panel.add(scene.component, BorderLayout.CENTER);
    // panel.add(buildSlider(), BorderLayout.SOUTH);
    
    frame.add(panel)
    frame.setSize(600, 600)
    frame.setVisible(true)
    frame
  }
}
