package kip.md

import java.awt.{Color, Frame, BorderLayout}
import javax.swing.JPanel
import java.awt.event._
import kip.graphics.{Scene, Bounds3d, GfxGL}
import kip.math.{Vec3, Quaternion}
import kip.util.Interpreter

case class Appearance(radius: Double, color: Color, sphereRes: Int = 2)
case class Particle(pt: Pt, appear: Appearance)


object Visualizer {
  def main(args: Array[String]) {
    val viz = new Visualizer()
    viz.setParticles(Seq(Particle(new Pt {var x=0.2; var y=0.8; var z=1.0}, Appearance(0.02, Color.RED))))
    Interpreter.start(("molviz", viz))
  }
}

class Visualizer {
  
  var particles = Seq[Particle]()
  var bds = Bounds3d(Vec3(0,0,0), Vec3(1,1,1))
  
  def setParticles(particles: Seq[Particle]) {
    this.particles = particles
  }

  def setBounds(bds: Bounds3d) {
    this.bds = bds;
  }

  val scene = new Scene() {
    def drawContent(gfx: GfxGL) {
      gfx.perspective3d(bds, rotation, translation*(bds.hi-bds.lo).norm)
      gfx.setColor(Color.GREEN)
      gfx.drawCuboid(bds)
      for (p <- particles) {
        gfx.setColor(p.appear.color)
        gfx.drawSphere(p.pt.toVec, p.appear.radius, p.appear.sphereRes)
      }
      gfx.ortho2dPixels()
      gfx.setColor(Color.RED)
      val pixOffset = 4
      gfx.rasterString("%d atoms".format(particles.size), pixOffset, pixOffset)
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
    frame.setSize(300, 300)
    frame.setVisible(true)
    frame
  }
}
