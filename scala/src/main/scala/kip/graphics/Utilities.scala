package kip.graphics

import javax.swing._
import java.awt.image.BufferedImage;
import javax.imageio.ImageIO;
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
  
  def writeImage(im: BufferedImage, file: String) {
    ImageIO.write(im, "png", new File(file));
  }
  
}
