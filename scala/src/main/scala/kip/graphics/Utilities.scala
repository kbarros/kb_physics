package kip.graphics

import javax.swing._
import java.awt.image.BufferedImage;
import javax.imageio.ImageIO;
import javax.media.opengl.{GL, GLCanvas}
import java.io.File


object Utilities {

  def captureJComponentImage(component: JComponent, width: Int, height: Int): BufferedImage = {
    val oldOpaque = component.isOpaque();
    val oldWidth = component.getWidth();
    val oldHeight = component.getHeight();
    component.setOpaque(true);
    component.setSize(width, height);
    
    val image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
    val g2d = image.createGraphics();
    component.paint(g2d);
    g2d.dispose();
    
    component.setOpaque(oldOpaque);
    component.setSize(oldWidth, oldHeight);
    return image;
  }
  

  // !! Use instead com.sun.opengl.util.Screenshot.readToBufferedImage() !!
  // -----------------------------------------------------------
  // The following two functions have been adapted from a posting by
  // keving on the JOGL board. The original post can be found here:
  // http://www.javagaming.org/cgi-bin/JGNetForums/YaBB.cgi?board=jogl;action=display;num=1091095319;start=2#2
  // -----------------------------------------------------------
  def copyFrame(canvas: GLCanvas): BufferedImage = { // copies the Frame to an integer array
    val w = canvas.getSize().width; // get the canvas' dimensions
    val h = canvas.getSize().height;
 
    val pixelsRGB = java.nio.ByteBuffer.allocateDirect(w * h * 3)
    val gl = canvas.getGL(); // acquire our GL Object
    
    // read the Frame back into our ByteBuffer
    gl.glReadBuffer(GL.GL_BACK);
    gl.glPixelStorei(GL.GL_PACK_ALIGNMENT, 1);
    gl.glReadPixels(0, 0, w, h, GL.GL_RGB, GL.GL_UNSIGNED_BYTE, pixelsRGB);
    
    val pixelInts = new Array[Int](w * h)
    
    // Convert RGB bytes to ARGB ints with no transparency. Flip image vertically by reading the
    // rows of pixels in the byte buffer in reverse - (0,0) is at bottom left in OpenGL.
    var p = w * h * 3; // Points to first byte (red) in each row.
    var q = 0   // Index into ByteBuffer
    var i = 0   // Index into target int[]
    val w3 = w*3;    // Number of bytes in each row
    for (row <- 0 until h) {
      p -= w3
      q = p
      for (col <- 0 until w) {
        val iR = pixelsRGB.get(q);
        q += 1
        val iG = pixelsRGB.get(q);
        q += 1
        val iB = pixelsRGB.get(q);
        q += 1
        pixelInts(i) = 0xFF000000 | ((iR & 0x000000FF) << 16) | ((iG & 0x000000FF) << 8) | (iB & 0x000000FF);
        i += 1
      }
    }
    val bufferedImage = new BufferedImage(w, h, BufferedImage.TYPE_INT_ARGB)
    bufferedImage.setRGB(0, 0, w, h, pixelInts, 0, w)
    bufferedImage
  }
}
