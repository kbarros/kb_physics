package kip.graphics

import java.awt.Frame
import java.awt.event.{WindowAdapter, WindowEvent}
import javax.media.opengl.{GL, GL2, GLAutoDrawable, GLEventListener}
import javax.media.opengl.fixedfunc.{GLMatrixFunc}
import javax.media.opengl.awt.GLCanvas


object HelloWorld {
  def main(args: Array[String]) {
    val frame = new Frame("JOGL HelloWorld")
    val canvas = new GLCanvas()
    canvas.addGLEventListener(new HelloWorld())
    frame.add(canvas)
    frame.setSize(300, 300)
    frame.setVisible(true)
    frame.addWindowListener(new WindowAdapter() {
      override def windowClosing(e: WindowEvent) {
	System.exit(0);
      }
    })
  }
}

class HelloWorld extends GLEventListener {
  def init(drawable: GLAutoDrawable) {
  }

  def reshape(drawable: GLAutoDrawable, x: Int, y: Int, width: Int, height: Int) {
  }

  def display(drawable: GLAutoDrawable) {
    val gl = drawable.getGL().getGL2()
    gl.glClearColor(0.0f, 0.0f, 0.0f, 0.0f)
    gl.glClear(GL.GL_COLOR_BUFFER_BIT);
    gl.glColor3f(1.0f, 1.0f, 1.0f)

    gl.glMatrixMode(GLMatrixFunc.GL_PROJECTION)
    gl.glLoadIdentity()
    gl.glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)
    gl.glMatrixMode(GLMatrixFunc.GL_MODELVIEW)
    gl.glLoadIdentity()
    
    gl.glBegin(GL2.GL_POLYGON)
    gl.glVertex2f(-0.5f, -0.5f)
    gl.glVertex2f(-0.5f, 0.5f)
    gl.glVertex2f(0.5f, 0.5f)
    gl.glVertex2f(0.5f, -0.5f)
    gl.glEnd()
    gl.glFlush()
  }

  def dispose(drawable: GLAutoDrawable) {
  }
}

