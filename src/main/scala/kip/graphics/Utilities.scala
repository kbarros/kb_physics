package kip.graphics

import javax.swing._
import java.awt.image.BufferedImage;
import javax.imageio.ImageIO;
import javax.media.opengl.{GL, GLCanvas}
import java.io.File


object Utilities {

  // Captures the image of a JComponent object.
  // Note that GLCanvas is not a subclass of JComponent; to capture the image of GLCanvas it is necessary
  // to use com.sun.opengl.util.Screenshot.readToBufferedImage instead. See Scene.display() for an example.
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


  // Use com.sun.opengl.util.Screenshot.readToBufferedImage() instead
  // See git history for definition.
  // 6/24/2011
  // def copyFrame(canvas: GLCanvas): BufferedImage = { ... }
}
