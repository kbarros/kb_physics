package kip.graphics

import java.awt.{Component, BorderLayout, Color}
import java.awt.event.{WindowAdapter, WindowEvent, MouseEvent}
import javax.swing.{JPanel, BorderFactory}
import javax.swing.event.{MouseInputAdapter, MouseInputListener}
import javax.media.opengl.{GL, GLAutoDrawable, GLEventListener, GLCapabilities}
import javax.media.opengl.fixedfunc.{GLLightingFunc}
import javax.media.opengl.awt.GLCanvas
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
  private var _capturingNow: Boolean = false
  private var _captureImage: java.awt.image.BufferedImage = _
  
  val (canvas, component) = initialize()
  
  def drawContent(gfx: GfxGL)
  
  def captureImage(): java.awt.image.BufferedImage = {
    _capturingNow = true
    _captureImage = null
    canvas.display()
    _captureImage
  }
  
  def display() {
    canvas.display()
  }
  
  def initialize(): (GLCanvas, JPanel) = {
    val capabilities = new GLCapabilities(null)
    capabilities.setSampleBuffers(true)
    capabilities.setNumSamples(4)
    val canvas = new GLCanvas(capabilities)

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
        val gl = drawable.getGL().getGL2()
        gl.glClearColor(1.0f, 1.0f, 1.0f, 0.0f)
        gl.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        gl.glEnable(GL.GL_BLEND)
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)
        gl.glEnable(GL.GL_LINE_SMOOTH)
        gl.glShadeModel(GLLightingFunc.GL_SMOOTH)
        drawContent(new GfxGL(drawable))
        gl.glFlush()
        
        if (_capturingNow) {
          _capturingNow = false
          _captureImage = com.jogamp.opengl.util.awt.Screenshot.readToBufferedImage(canvas.getWidth(), canvas.getHeight())
        }
      }
      
      def dispose(drawable: GLAutoDrawable) {
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
