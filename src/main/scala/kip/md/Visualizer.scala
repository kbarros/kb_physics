package kip.md

import java.awt.{Color, Frame, BorderLayout}
import javax.swing.JPanel
import java.awt.event._
import kip.graphics.{Scene, Bounds3d, GfxGL}
import kip.math.{Vec3, Quaternion}
import kip.util.Interpreter

object Visualizer {
  case class Sphere(pt: Vec3, radius: Double, color: Color, resolution: Int = 2)
  case class Wall(norm: Vec3, pos: Vec3, color: Color)
  case class Arrow(from: Vec3, to: Vec3, normal: Vec3, color: Color)
  case class Path(pts: Seq[Vec3], color: Color)
}

class Visualizer(sizew: Int = 600, sizeh: Int = 600) {
  
  var particles = Seq[Visualizer.Sphere]()
  var walls = Seq[Visualizer.Wall]()
  var arrows = Seq[Visualizer.Arrow]()
  var paths = Seq[Visualizer.Path]()
  var bds = Bounds3d(Vec3(0,0,0), Vec3(1,1,1))
  var rasterString = ""
  var arrowHeadSize = 1.0
  
  def setBounds(bds: Bounds3d) {
    this.bds = bds
  }
  
  def setParticles(particles: Seq[Visualizer.Sphere]) {
    this.particles = particles
  }
  
  def setWalls(walls: Seq[Visualizer.Wall]) {
    this.walls = walls
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
      for (w <- walls) {
        gfx.setColor(w.color)
        // TODO: construct proper boundary for wall, perhaps in Volume?
        // c.f. Grid3DSliceView.draw()
        val len = (bds.hi-bds.lo).norm // typical length of volume
        val tangent = Vec3(w.norm.y, -w.norm.x, 0) // project norm onto (x,y) and rotate by Pi/2
        val p1 = w.pos - tangent*len
        val p2 = w.pos + tangent*len
        gfx.drawLines(p1, p2)
      }
      
      // gfx.setMultisampling(true)
      for (s <- arrows) {
        gfx.setColor(s.color)
        val x = (s.to - s.from).normalize
        val t = s.normal cross x         
        val mid = s.to - x*arrowHeadSize 
        
        val arrowHeadWidth = arrowHeadSize / 2
        val arrowShaftWidth = arrowHeadSize / 8
        
        if ((s.to - s.from).norm > arrowHeadSize) {
          val w = t*arrowShaftWidth
          gfx.drawQuads(s.from+w, s.from-w, mid-w, mid+w)
        }
        
        val v0 = s.to
        val v1 = mid + t*(arrowHeadWidth)
        val v2 = mid - t*(arrowHeadWidth)
        gfx.drawTriangles(v0, v1, v2)
      }
      
      gfx.gl.glLineWidth(2)
      for (p <- paths) {
        gfx.setColor(p.color)
        gfx.drawLineStrip(p.pts: _*)
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
    frame.setSize(sizew, sizeh)
    frame.setVisible(true)
    frame
  }
}
