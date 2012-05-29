package kip.graphics

import java.awt.Frame
import java.awt.event.{WindowAdapter, WindowEvent}
import javax.media.opengl.{GL, GL2, GLAutoDrawable, GLEventListener}
import javax.media.opengl.fixedfunc.{GLMatrixFunc}
import javax.media.opengl.awt.GLCanvas
import java.nio.ByteBuffer;
import javax.swing.SwingUtilities


object GridView {
  class Data(val w: Int, val h: Int, val color: Int => java.awt.Color)
  class ArrayData(w: Int, h: Int, a: Array[Double], cg: ColorGradient) extends Data(w, h, i => cg.interpolate(a(i))) {
    if (w*h != a.size) println("Warning: array.size (%d) != w*h (%dx%d)".format(a.size, w, h))
  }
  
  def main(args: Array[String]) {
    val gridView = new GridView()
    val data = {
      val w = 256
      val data = Array.tabulate(w, w) { (i, j) => (i ^ j).toDouble }
      new GridView.ArrayData(w, w, data.flatten, ColorGradient.blueRed(lo=0, hi=256))
    }
    gridView.display(data)
    
    Utilities.frame(gridView.canvas, w=300, h=300, title="GridView Test")
  }
}


class GridView {
  var data: GridView.Data = new GridView.Data(1, 1, _ => java.awt.Color.BLACK)
  
  val listener: GLEventListener = new GLEventListener {
    def init(drawable: GLAutoDrawable) {
    }

    def reshape(drawable: GLAutoDrawable, x: Int, y: Int, width: Int, height: Int) {
    }

    def display(drawable: GLAutoDrawable) {
      val gl = drawable.getGL().getGL2()
      
      gl.glClearColor(1f, 1f, 1f, 0f)
      gl.glClear(GL.GL_COLOR_BUFFER_BIT);
      
      gl.glMatrixMode(GLMatrixFunc.GL_PROJECTION)
      gl.glLoadIdentity()
      gl.glOrtho(0, 1, 0, 1, -1, 1)
      gl.glMatrixMode(GLMatrixFunc.GL_MODELVIEW)
      gl.glLoadIdentity()
      
      // allocate and bind texture
      gl.glEnable(GL.GL_TEXTURE_2D)
      val textures = new Array[Int](1)
      gl.glGenTextures(textures.size, textures, 0)
      gl.glBindTexture(GL.GL_TEXTURE_2D, textures(0))
      gl.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_NEAREST /* GL.GL_LINEAR */)
      gl.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_NEAREST /* GL.GL_LINEAR */)
      gl.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_S, GL.GL_CLAMP_TO_EDGE)
      gl.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_T, GL.GL_CLAMP_TO_EDGE)
      
      // load texture
      val buffer = com.jogamp.common.nio.Buffers.newDirectByteBuffer(4*data.w*data.h)
      buffer.clear()
      for (i <- 0 until data.w*data.h) {
        val c = data.color(i)
        buffer.put(c.getRed().toByte)
        buffer.put(c.getGreen().toByte)
        buffer.put(c.getBlue().toByte)
        buffer.put(c.getAlpha().toByte)
      }
      buffer.flip()
      gl.glTexImage2D(
        GL.GL_TEXTURE_2D, 0, GL.GL_RGBA, data.w, data.h,
        0, GL.GL_RGBA, GL.GL_UNSIGNED_BYTE, buffer)
      
      // display
      val lo = 0.02f
      val hi = 0.98f
      gl.glColor3f(1, 1, 1)
      gl.glBegin(GL2.GL_POLYGON)
      gl.glTexCoord2d(0, 0)
      gl.glVertex2f(lo, lo)
      gl.glTexCoord2d(0, 1)
      gl.glVertex2f(lo, hi)
      gl.glTexCoord2d(1, 1)
      gl.glVertex2f(hi, hi)
      gl.glTexCoord2d(1, 0)
      gl.glVertex2f(hi, lo)
      gl.glEnd()
      gl.glColor3f(0, 0, 0)
      gl.glBegin(GL.GL_LINE_LOOP)
      gl.glVertex2f(lo, lo)
      gl.glVertex2f(lo, hi)
      gl.glVertex2f(hi, hi)
      gl.glVertex2f(hi, lo)
      gl.glEnd()
      
      // clean up
      gl.glDeleteTextures(textures.size, textures, 0);
      gl.glDisable(GL.GL_TEXTURE_2D);		
      
      gl.glFlush()
    }
    
    def dispose(drawable: GLAutoDrawable) {
    }
  }
  
  val canvas: GLCanvas = {
    val canvas = new GLCanvas()
    canvas.addGLEventListener(listener)
    canvas
  }
  
  def display(d: GridView.Data) {
    data = d
    Utilities.swingInvokeLater(canvas.display)
  }
  
  def clear() {
    display(new GridView.Data(1, 1, _ => java.awt.Color.BLACK))
  }
}

