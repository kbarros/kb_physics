package kip.graphics


import java.awt.{Color, Frame, BorderLayout}
import javax.swing.JPanel
import java.awt.event._
import kip.math.{Vec3, Quaternion}


object RetainedScene {
  trait Drawable {
    def draw(gfx: GfxGL)
  }
  class Sphere(var pt: Vec3, var radius: Double, var color: Color = Color.BLUE, var resolution: Int = 2) extends Drawable {
    def draw(gfx: GfxGL) {
      gfx.setColor(color)
      gfx.drawSphere(pt, radius, resolution)
    }
  }
  
  class Arrow(var origin: Vec3, var delta: Vec3, var width: Double,
              var color1: Color = Color.BLUE, var color2: Color = Color.MAGENTA) extends Drawable {
    def draw(gfx: GfxGL) {
      val a0 = delta.normalize 
      val a1 = {
        val alignedz = a0.z > (1 - 1e-6) 
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
      
      gfx.setColor(color1)
      // bottom cap
      gfx.drawQuads(x0 + (+a1+a2)*w0, x0 + (+a1-a2)*w0, x1 + (-a1-a2)*w1, x1 + (-a1+a2)*w1)      
      // base
      gfx.drawQuads(x0 + (+a1-a2)*w0, x0 + (+a1+a2)*w0, x1 + (+a1+a2)*w1, x1 + (+a1-a2)*w1)
      gfx.drawQuads(x0 + (+a1+a2)*w0, x0 + (-a1+a2)*w0, x1 + (-a1+a2)*w1, x1 + (+a1+a2)*w1)
      gfx.drawQuads(x0 + (-a1+a2)*w0, x0 + (-a1-a2)*w0, x1 + (-a1-a2)*w1, x1 + (-a1+a2)*w1)
      gfx.drawQuads(x0 + (-a1-a2)*w0, x0 + (+a1-a2)*w0, x1 + (+a1-a2)*w1, x1 + (-a1-a2)*w1)
      // head
      gfx.setColor(color2)
      gfx.drawTriangles(x1 + (+a1-a2)*w1, x1 + (+a1+a2)*w1, x2)
      gfx.drawTriangles(x1 + (+a1+a2)*w1, x1 + (-a1+a2)*w1, x2)
      gfx.drawTriangles(x1 + (-a1+a2)*w1, x1 + (-a1-a2)*w1, x2)
      gfx.drawTriangles(x1 + (-a1-a2)*w1, x1 + (+a1-a2)*w1, x2)

      //      gfx.drawQuads(x0 + (+a1+a2)*w0, x0 + (-a1+a2)*w0, x1 + (-a1+a2)*w1, x1 + (+a1+a2)*w1)
//      gfx.drawQuads(x0 + (-a1+a2)*w0, x0 + (-a1-a2)*w0, x1 + (-a1-a2)*w1, x1 + (-a1+a2)*w1)
//      gfx.drawQuads(x0 + (-a1-a2)*w0, x0 + (+a1-a2)*w0, x1 + (+a1-a2)*w1, x1 + (-a1-a2)*w1)
    }
  }
  
  class Cuboid(bds: Bounds3d, color: Color = Color.GREEN) extends Drawable {
    def draw(gfx: GfxGL) {
      gfx.setColor(color)
      gfx.drawCuboid(bds)
    }
  }
}

class RetainedScene(var bds: Bounds3d, sizew: Int = 600, sizeh: Int = 600) {
  var drawables = Vector[RetainedScene.Drawable]()
  var rasterString = ""
  var arrowHeadSize = 1.0
  
  def display() {
    scene.display()
  }
  
  val scene = new Scene() {
    def drawContent(gfx: GfxGL) {
      gfx.perspective3d(bds, rotation, translation*(bds.hi-bds.lo).norm)
      
      // gfx.setMultisampling(true)
      for (d <- drawables) {
        d.draw(gfx)
      }
      
      gfx.ortho2dPixels()
      gfx.setColor(Color.RED)
      val pixOffset = 4
      gfx.rasterString(rasterString, pixOffset, pixOffset)
    }
  }
  
  val frame = buildFrame()
  
  def buildFrame(): Frame = {
    val frame = new Frame("Retained Scene")
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
