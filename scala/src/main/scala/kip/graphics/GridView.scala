package kip.graphics

import java.awt.Frame
import java.awt.event.{WindowAdapter, WindowEvent}
import javax.media.opengl.{GL, GLAutoDrawable, GLCanvas, GLEventListener}
import com.sun.opengl.util.BufferUtil;
import java.nio.ByteBuffer;


object GridView {
  case class Data(w: Int, h: Int, a: Array[Double]) {
    assert(w*h == a.size)
  }
  
  def main(args: Array[String]) {
    val gridView = new GridView()
    gridView.data = {
      val w = 256
      val data = Array.tabulate(w, w) { (i, j) => (i ^ j) / w.toDouble }
      GridView.Data(w, w, data.flatten)
    }
    
    val frame = new Frame("GridView Test")
    frame.add(gridView.canvas)
    frame.setSize(300, 300)
    frame.setVisible(true)
    frame.addWindowListener(new WindowAdapter() {
      override def windowClosing(e: WindowEvent) {
	System.exit(0);
      }
    })
  }
}


class GridView {
  var data: GridView.Data = GridView.Data(1, 1, Array(0d))
  
  val listener: GLEventListener = new GLEventListener {
    def init(drawable: GLAutoDrawable) {
    }

    def reshape(drawable: GLAutoDrawable, x: Int, y: Int, width: Int, height: Int) {
    }

    def display(drawable: GLAutoDrawable) {
      val gl = drawable.getGL()
      
      gl.glClearColor(1f, 1f, 1f, 0f)
      gl.glClear(GL.GL_COLOR_BUFFER_BIT);

      gl.glMatrixMode(GL.GL_PROJECTION)
      gl.glLoadIdentity()
      gl.glOrtho(0, 1, 0, 1, -1, 1)
      gl.glMatrixMode(GL.GL_MODELVIEW)
      gl.glLoadIdentity()
      
      gl.glEnable(GL.GL_TEXTURE_2D)
      
      // load texture
      val textures = new Array[Int](1)
      gl.glGenTextures(textures.size, textures, 0);
      gl.glBindTexture(GL.GL_TEXTURE_2D, textures(0))
      
      val buffer = BufferUtil.newByteBuffer(4*data.a.size);
      buffer.clear()
      for (i <- 0 until data.a.size) {
        val c = ColorGradient.interpolate(data.a(i))
        buffer.put(c.getRed().toByte)
        buffer.put(c.getGreen().toByte)
        buffer.put(c.getBlue().toByte)
        buffer.put(c.getAlpha().toByte)
      }
      buffer.flip()
      
      // TODO: move this init stuff up
      gl.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_LINEAR);
      gl.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_LINEAR);
      gl.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_S, GL.GL_CLAMP_TO_EDGE)
      gl.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_T, GL.GL_CLAMP_TO_EDGE)
      gl.glTexImage2D(
        GL.GL_TEXTURE_2D, 0, GL.GL_RGBA, data.w, data.h,
        0, GL.GL_RGBA, GL.GL_UNSIGNED_BYTE, buffer)
      
      // display
      val lo = 0.02f
      val hi = 0.98f
      gl.glColor3f(1, 1, 1)
      gl.glBegin(GL.GL_POLYGON)
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
    
    def displayChanged(drawable: GLAutoDrawable, modeChanged: Boolean, deviceChanged: Boolean) {
    }
  }

  val canvas: GLCanvas = {
    val canvas = new GLCanvas()
    canvas.addGLEventListener(listener)
    canvas
  }
  
}

