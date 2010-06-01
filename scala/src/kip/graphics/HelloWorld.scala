package kip.graphics

import java.awt.Frame
import java.awt.event.{WindowAdapter, WindowEvent}
import javax.media.opengl.{GL, GLAutoDrawable, GLCanvas, GLEventListener}


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
    val gl = drawable.getGL()
    gl.glClearColor(0.0f, 0.0f, 0.0f, 0.0f)
    gl.glClear(GL.GL_COLOR_BUFFER_BIT);
    gl.glColor3f(1.0f, 1.0f, 1.0f)
    gl.glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0)
    gl.glBegin(GL.GL_POLYGON)
    gl.glVertex2f(-0.5f, -0.5f)
    gl.glVertex2f(-0.5f, 0.5f)
    gl.glVertex2f(0.5f, 0.5f)
    gl.glVertex2f(0.5f, -0.5f)
    gl.glEnd()
    gl.glFlush()
  }

  def displayChanged(drawable: GLAutoDrawable, modeChanged: Boolean, deviceChanged: Boolean) {
  }
}

