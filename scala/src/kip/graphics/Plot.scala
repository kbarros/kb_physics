package kip.graphics

import java.awt.{BorderLayout, Frame, Color}
import java.awt.event.{WindowAdapter, WindowEvent, MouseEvent}
import javax.swing.{JPanel, BorderFactory}
import javax.swing.event.{MouseInputAdapter, MouseInputListener}
import javax.media.opengl.{GL, GLAutoDrawable, GLCanvas, GLEventListener}
import kip.math.{Vec3, Quaternion}


trait DragHandler extends MouseInputAdapter {
  var lastDrag: Option[java.awt.Point] = None
  
  override def mousePressed(e: MouseEvent) {
    if (!e.isPopupTrigger())
      lastDrag = Some(e.getPoint())
  }
  
  override def mouseReleased(e: MouseEvent) {
    lastDrag = None
  }
  
  override def mouseDragged(e: MouseEvent) {
    lastDrag foreach{ ld =>
      val dx = e.getX() - ld.x
      val dy = e.getY() - ld.y
      lastDrag = Some(e.getPoint())
      mouseDraggedDelta(dx, dy, e)
    }
  }
  
  def mouseDraggedDelta(dx: Int, dy: Int, e: MouseEvent)
}

/**
 *
 */
class Scene {
  val canvas = new GLCanvas()
  var rotation: Quaternion = Quaternion.fromAxisAngle(0, 0, 0)
  
  val mouse = new DragHandler {
    def mouseDraggedDelta(dx: Int, dy: Int, e: MouseEvent) {
      val rpp = 0.01 // radians per pixel
      val q = Quaternion.fromAxisAngle(dy*rpp, dx*rpp, 0)
      rotation = (q * rotation).normalize
      canvas.repaint()
    }
  }
  canvas.addMouseListener(mouse)
  canvas.addMouseMotionListener(mouse)
  
  canvas.addGLEventListener(new GLEventListener {
    def init(drawable: GLAutoDrawable) {
    }

    def reshape(drawable: GLAutoDrawable, x: Int, y: Int, width: Int, height: Int) {
    }

    def display(drawable: GLAutoDrawable) {
      val gl = drawable.getGL()
      gl.glClearColor(1.0f, 1.0f, 1.0f, 0.0f)
      gl.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
      gl.glEnable(GL.GL_BLEND)
      gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)
      gl.glEnable(GL.GL_LINE_SMOOTH)
      gl.glShadeModel(GL.GL_SMOOTH)
      
      val bds = Bounds3d(Vec3(0,0,0), Vec3(100, 100, 100))
      val gfx = new GfxGL(drawable)
//      gfx.ortho2d(bds)
      gfx.perspective3d(bds, rotation)
      gfx.setColor(Color.BLUE)
      gfx.drawLines(Vec3(0,0,0), Vec3(50,80,50))
      gfx.setColor(Color.GREEN)
      gfx.drawCuboid(bds)
      gfx.drawSphere(Vec3(20, 20, 20), 10)
      
      gl.glFlush()
    }

    def displayChanged(drawable: GLAutoDrawable, modeChanged: Boolean, deviceChanged: Boolean) {
    }
  })
    

  val component = new JPanel(new BorderLayout())
  component.setBorder(BorderFactory.createCompoundBorder(
    BorderFactory.createEmptyBorder(4, 4, 4, 4),
    BorderFactory.createLineBorder(Color.GRAY)))
  component.add(canvas)

}


object Plot {
  def main(args: Array[String]) {
    val scene = new Scene()
    
    val frame = new Frame("JOGL HelloWorld2")
    frame.add(scene.component)
    frame.setSize(300, 300)
    frame.setVisible(true)
    frame.addWindowListener(new WindowAdapter() {
      override def windowClosing(e: WindowEvent) {
	System.exit(0)
      }
    })
  }
}

