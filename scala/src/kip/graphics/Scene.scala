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
abstract class Scene {
  var rotation: Quaternion = Quaternion.fromAxisAngle(0, 0, 0)
  var translation: Vec3 = Vec3(0, 0, 0)
  
  val (canvas, component) = initialize()
  
  def drawContent(gfx: GfxGL)
  
  // def invokeInSwing(f: => Unit) {
  //   java.awt.EventQueue.invokeAndWait(new Runnable {
  //     def run() = f
  //   })
  // }
  
  def triggerRepaint() {
    canvas.repaint()
  }
  
  def initialize(): (GLCanvas, JPanel) = {
    val canvas = new GLCanvas()
    val mouse = new DragHandler {
      def mouseDraggedDelta(dx: Int, dy: Int, e: MouseEvent) {
        val dpp = 0.003 // dimensionless displacement per pixel
        val rpp = 0.01  // radians per pixel
        if (e.isShiftDown) {
          translation += Vec3(dx, -dy, 0) * dpp
        }
        else if (e.isMetaDown) {
          translation += Vec3(0, 0, -dy)*dpp
          // val q = Quaternion.fromAxisAngle(0, 0, -dx*rpp)
          // rotation = (q * rotation).normalize
        }
        else {
          val q = Quaternion.fromAxisAngle(dy*rpp, dx*rpp, 0)
          rotation = (q * rotation).normalize
        }
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
        drawContent(new GfxGL(drawable))
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
    
    (canvas, component)
  }

}

