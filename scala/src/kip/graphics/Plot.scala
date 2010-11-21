package kip.graphics

import kip.math.{Vec3}
import java.awt.{Color, Frame}
import java.awt.event.{WindowAdapter, WindowEvent}


object Plot {
  def main(args: Array[String]) {
    val scene = new Scene() {
      def drawContent(gfx: GfxGL) {
        val bds = Bounds3d(Vec3(0,0,0), Vec3(100, 100, 100))
        gfx.perspective3d(bds, rotation, translation)
        gfx.setColor(Color.BLUE)
        gfx.drawLines(Vec3(0,0,0), Vec3(50,80,50))
        gfx.setColor(Color.GREEN)
        gfx.drawCuboid(bds)
        gfx.drawSphere(Vec3(20, 20, 20), 10, 0)
        gfx.ortho2d(Bounds3d(Vec3(0,0,0), Vec3(1,1,1)))
        gfx.setColor(Color.RED)
        gfx.rasterString("hello worl", 0.1, 0.1)
      }
    }
    
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

