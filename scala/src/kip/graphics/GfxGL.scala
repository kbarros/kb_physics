package kip.graphics

import java.awt.Color
import javax.media.opengl.{GL, GLAutoDrawable}
import javax.media.opengl.glu.GLU
import com.sun.opengl.util.GLUT

import kip.math.Math._
import kip.math.{Vec3, Quaternion}


case class Bounds3d(lo: Vec3, hi: Vec3) {
  def width  = hi.x - lo.x
  def height = hi.y - lo.y
  def depth  = hi.z - lo.z
  def center = (lo+hi) / 2
}

class GfxGL(glDrawable: GLAutoDrawable) {
  val gl   = glDrawable.getGL()
  val glu  = new GLU()
  val gluq = glu.gluNewQuadric();
  val glut = new GLUT()
  val pixWidth  = glDrawable.getWidth()
  val pixHeight = glDrawable.getHeight()
  
  val FONT = GLUT.BITMAP_8_BY_13
  val FONT_HEIGHT = 13 // pixels

  
  def ortho2d(bds: Bounds3d) {
    gl.glDisable(GL.GL_DEPTH_TEST)
    gl.glDisable(GL.GL_COLOR_MATERIAL)
    gl.glDisable(GL.GL_LIGHTING)
    
    // set the projection & modelview matrices
    gl.glMatrixMode(GL.GL_PROJECTION)
    gl.glLoadIdentity()
    glu.gluOrtho2D(bds.lo.x, bds.hi.x, bds.lo.y, bds.hi.y)
    gl.glMatrixMode(GL.GL_MODELVIEW)
    gl.glLoadIdentity()
  }

  def perspective3d(bds: Bounds3d, rotation: Quaternion) {
    gl.glEnable(GL.GL_DEPTH_TEST)
    gl.glEnable(GL.GL_COLOR_MATERIAL)
    gl.glEnable(GL.GL_LIGHTING)
    gl.glLightModeli(GL.GL_LIGHT_MODEL_TWO_SIDE, GL.GL_TRUE)
    
    // get the corner to corner distance of the view bounds cuboid
    val len = hypot(bds.width, bds.height, bds.depth)
    
    // set the projection matrix
    gl.glMatrixMode(GL.GL_PROJECTION)
    gl.glLoadIdentity()
    val fovY = 35
    val aspect = pixWidth / pixHeight.toDouble
    val zNear = 0.1*len
    val zFar = 10*len
    glu.gluPerspective(fovY, aspect, zNear, zFar)
    
    // set the modelview matrix
    gl.glMatrixMode(GL.GL_MODELVIEW)
    gl.glLoadIdentity()
    initializeLights()
    // each sequential operation multiplies the modelview transformation matrix
    // from the left. operations on the scene object occur in reverse order from
    // their specification.
    // step (3): move object away from camera
    gl.glTranslated(0, 0, -1.5*len)
    // step (2): rotate object about zero
    gl.glMultMatrixd(rotation.toGLRotationMatrix, 0)
    // step (1): move object to its center
    val center = bds.center
    gl.glTranslated(-center.x, -center.y, -center.z)
  }

  def setLineSmoothing(b: Boolean) {
    if (b)
      gl.glEnable(GL.GL_LINE_SMOOTH)
    else
      gl.glDisable(GL.GL_LINE_SMOOTH)
  }
  
  def setColor(color: Color) {
    gl.glColor4fv(color.getComponents(null), 0)
  }
  
  def drawPoints(ps: Vec3*) {
    gl.glBegin(GL.GL_POINTS)
    for (p <- ps) {
      gl.glVertex3d(p.x, p.y, p.z)
    }
    gl.glEnd()
  }
  
  def drawLines(ps: Vec3*) {
    gl.glBegin(GL.GL_LINES)
    for (p <- ps) {
      gl.glVertex3d(p.x, p.y, p.z)
    }
    gl.glEnd();
  }
  
  def drawCuboid(bds: Bounds3d) {
    gl.glDisable(GL.GL_LIGHTING)
    gl.glBegin(GL.GL_LINES)
    val xs = Array(bds.lo.x, bds.hi.x)
    val ys = Array(bds.lo.y, bds.hi.y)
    val zs = Array(bds.lo.z, bds.hi.z)
    for (i <- 0 until 2;
         j <- 0 until 2;
         k <- 0 until 2;
         if (i + j + k) % 2 == 0) {
	    gl.glVertex3d(xs(i),   ys(j),   zs(k))
	    gl.glVertex3d(xs(1-i), ys(j),   zs(k))
	    gl.glVertex3d(xs(i),   ys(j),   zs(k))
	    gl.glVertex3d(xs(i),   ys(1-j), zs(k))
	    gl.glVertex3d(xs(i),   ys(j),   zs(k))
	    gl.glVertex3d(xs(i),   ys(j),   zs(1-k))
    }
    gl.glEnd()
    gl.glEnable(GL.GL_LIGHTING)
  }
  
  def drawSphere(center: Vec3, radius: Double) {
    gl.glPushMatrix()
    gl.glTranslated(center.x, center.y, center.z)
    glu.gluSphere(gluq, radius, 8, 8)
    gl.glPopMatrix()
  }
  
  def stringWidth(str: String) {
    glut.glutBitmapLength(FONT, str)
  }
	
  def stringHeight(str: String) = FONT_HEIGHT
  
  def rasterString(str: String, x: Double, y: Double) {
    gl.glPushMatrix()
    gl.glRasterPos2d(x, y)
    glut.glutBitmapString(FONT, str) 
    gl.glPopMatrix()
  }
  
  def initializeLights() {
    gl.glEnable(GL.GL_LIGHT1)
    gl.glLightfv(GL.GL_LIGHT1, GL.GL_AMBIENT,  Array(0.2f,0.2f,0.2f,0.2f), 0)
    gl.glLightfv(GL.GL_LIGHT1, GL.GL_DIFFUSE,  Array(0.9f,0.9f,0.9f,0.9f), 0)
    gl.glLightfv(GL.GL_LIGHT1, GL.GL_POSITION, Array(1,0.5f,1,0), 0)
  }
}
