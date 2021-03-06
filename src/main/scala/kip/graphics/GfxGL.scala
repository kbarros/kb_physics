package kip.graphics

import scala.collection.mutable.LinkedList
import java.awt.Color
import javax.media.opengl.{GL, GL2ES1, GL2, GLAutoDrawable}
import javax.media.opengl.fixedfunc.{GLLightingFunc, GLMatrixFunc}
import javax.media.opengl.glu.GLU
import com.jogamp.opengl.util.gl2.GLUT

import kip.math.Math._
import kip.math.{Vec3, Quaternion}



case class Bounds3d(lo: Vec3, hi: Vec3) {
//  def dx = hi.x - lo.x
//  def dy = hi.y - lo.y
//  def dz = hi.z - lo.z
  def center = (lo+hi) / 2
}

class GfxGL(glDrawable: GLAutoDrawable) {
  import GL._
  import GLLightingFunc._
  import GLMatrixFunc._
  import GL2ES1._
  import GL2._ // This should include previous imports, but doesn't.
  
  val gl   = glDrawable.getGL().getGL2()
  val glu  = new GLU()
  val gluq = glu.gluNewQuadric()
  val glut = new GLUT()
  val pixWidth  = glDrawable.getWidth()
  val pixHeight = glDrawable.getHeight()
  
  val FONT = GLUT.BITMAP_8_BY_13
  val FONT_HEIGHT = 13 // pixels

  
  def ortho2dPixels() {
    ortho2d(Bounds3d(Vec3(0,0,0), Vec3(pixWidth,pixHeight,0)))
  }
  
  def ortho2d(bds: Bounds3d) {
    gl.glDisable(GL_DEPTH_TEST)
    gl.glDisable(GL_COLOR_MATERIAL)
    gl.glDisable(GL_LIGHTING)
    
    // set the projection & modelview matrices
    gl.glMatrixMode(GL_PROJECTION)
    gl.glLoadIdentity()
    glu.gluOrtho2D(bds.lo.x, bds.hi.x, bds.lo.y, bds.hi.y)
    gl.glMatrixMode(GL_MODELVIEW)
    gl.glLoadIdentity()
  }

  def perspective3d(bds: Bounds3d, rotation: Quaternion, translation: Vec3, cameraDistance: Double = 1.5) {
    gl.glEnable(GL_DEPTH_TEST)
    gl.glEnable(GL_COLOR_MATERIAL)
    gl.glEnable(GL_LIGHTING)
    gl.glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE)
    
    // get the corner to corner distance of the view bounds cuboid
    val len = (bds.hi - bds.lo).norm
    
    // set the projection matrix
    gl.glMatrixMode(GL_PROJECTION)
    gl.glLoadIdentity()
    val fovY = 35
    val aspect = pixWidth / pixHeight.toDouble
    val zNear = 0.1*len
    val zFar = 10*len
    glu.gluPerspective(fovY, aspect, zNear, zFar)
    
    // set the modelview matrix
    gl.glMatrixMode(GL_MODELVIEW)
    gl.glLoadIdentity()
    initializeLights()
    // each sequential operation multiplies the modelview transformation matrix
    // from the left. operations on the scene object occur in reverse order from
    // their specification.
    // step (4): translate object according to user
    gl.glTranslated(translation.x, translation.y, translation.z)
    // step (3): move object away from camera
    gl.glTranslated(0, 0, -cameraDistance*len)
    // step (2): rotate object about zero
    gl.glMultMatrixd(rotation.toGLRotationMatrix, 0)
    // step (1): move object to its center
    val center = bds.center
    gl.glTranslated(-center.x, -center.y, -center.z)
  }

  def setLineSmoothing(b: Boolean) {
    if (b)
      gl.glEnable(GL_LINE_SMOOTH)
    else
      gl.glDisable(GL_LINE_SMOOTH)
  }
  
  def setMultisampling(b: Boolean) {
    if (b)
      gl.glEnable(GL_MULTISAMPLE)
    else
      gl.glDisable(GL_MULTISAMPLE)
  }
  
  def setSmoothShading(b: Boolean) {
    if (b)
      gl.glShadeModel(GL_SMOOTH)
    else
      gl.glShadeModel(GL_FLAT)
  }
  
  def setColor(color: Color) {
    gl.glColor4fv(color.getComponents(null), 0)
  }
  
  def drawPoints(ps: Vec3*) {
    gl.glDisable(GL_LIGHTING)
    gl.glBegin(GL_POINTS)
    for (p <- ps) {
      gl.glVertex3d(p.x, p.y, p.z)
    }
    gl.glEnd()
    gl.glEnable(GL_LIGHTING)
  }
  
  def drawLines(ps: Vec3*) {
    gl.glDisable(GL_LIGHTING)
    gl.glBegin(GL_LINES)
    for (p <- ps) {
      gl.glVertex3d(p.x, p.y, p.z)
    }
    gl.glEnd();
    gl.glEnable(GL_LIGHTING)
  }
  
  def drawLineStrip(ps: Vec3*) {
    gl.glDisable(GL_LIGHTING)
    gl.glBegin(GL_LINE_STRIP)
    for (p <- ps) {
      gl.glVertex3d(p.x, p.y, p.z)
    }
    gl.glEnd();
    gl.glEnable(GL_LIGHTING)
  }
      
  def drawLineStrip(ps: IndexedSeq[Vec3], colors: IndexedSeq[Color]) {
    require(colors.size == ps.size-1)
    gl.glDisable(GL_LIGHTING)
    gl.glBegin(GL_LINE_STRIP)
    gl.glVertex3d(ps(0).x, ps(0).y, ps(0).z)
    for (i <- 0 until colors.size) {
      setColor(colors(i))
      gl.glVertex3d(ps(i+1).x, ps(i+1).y, ps(i+1).z)
    }
    gl.glEnd()
    gl.glEnable(GL_LIGHTING)
  }

  def drawTriangles(ps: Vec3*) {
    gl.glBegin(GL_TRIANGLES)
    for (qs <- ps.grouped(3)) {
      val n = ((qs(1) - qs(0)) cross (qs(2) - qs(0))).normalize
      gl.glNormal3d(n.x, n.y, n.z)
      for (q <- qs) {
        gl.glVertex3d(q.x, q.y, q.z)
      }
    }
    gl.glEnd();
  }

  def drawQuads(ps: Vec3*) {
    val ps2 = ps.grouped(4).flatMap {qs =>
      Seq(qs(0), qs(1), qs(2),
          qs(0), qs(2), qs(3))
    }
    drawTriangles(ps2.toSeq:_*)
  }

  def drawTriangleStrip(ps: IndexedSeq[Vec3], colors: IndexedSeq[Color]) {
    require(colors.size == ps.size-2)
    gl.glDisable(GL_LIGHTING)
    gl.glBegin(GL_TRIANGLE_STRIP)
    gl.glNormal3d(0, 0, 1)
    gl.glVertex3d(ps(0).x, ps(0).y, ps(0).z)
    gl.glVertex3d(ps(1).x, ps(1).y, ps(1).z)
    for (i <- 2 until ps.size) {
      setColor(colors(i-2))
      gl.glVertex3d(ps(i).x, ps(i).y, ps(i).z)
    }
    gl.glEnd()
  }
  
  def drawCuboid(bds: Bounds3d) {
    gl.glDisable(GL_LIGHTING)
    gl.glBegin(GL_LINES)
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
    gl.glEnable(GL_LIGHTING)
  }
  
  
  case class Triangle(v1: Vec3, v2: Vec3, v3: Vec3)
  val octahedron = {
    val a1 = Vec3(0, 0, 1)
    val b1 = Vec3(1, 0, 0)
    val b2 = Vec3(0, 1, 0)
    val b3 = Vec3(-1, 0, 0)
    val b4 = Vec3(0, -1, 0)
    val c1 = Vec3(0, 0, -1)
    LinkedList(
      Triangle(a1, b1, b2),
      Triangle(a1, b2, b3),
      Triangle(a1, b3, b4),
      Triangle(a1, b4, b1),
      Triangle(c1, b2, b1),
      Triangle(c1, b3, b2),
      Triangle(c1, b4, b3),
      Triangle(c1, b1, b4)
    )
  }
  
  val sphereCache = new Array[Seq[Triangle]](5)
  sphereCache(0) = octahedron
  
  def drawSphere(center: Vec3, radius: Double, subdivisions: Int) {
    if (subdivisions >= sphereCache.size)
      throw new IllegalArgumentException("Cannot subdivide sphere with n=%d".format(subdivisions))
    
    if (sphereCache(subdivisions) == null) {
      def subdivide(t: Triangle, n: Int): LinkedList[Triangle] = {
        if (n == 0)
          LinkedList(t)
        else {
          val v12 = (t.v1 + t.v2).normalize
          val v23 = (t.v2 + t.v3).normalize
          val v31 = (t.v3 + t.v1).normalize
          //     v23
          //v3 ________ v2
          //   \  /\  /
          // v31\/__\/ v12
          //     \  /
          //      \/ v1
          //       
          val triangles = LinkedList(
            Triangle(t.v1, v12, v31),
            Triangle(v12, t.v2, v23),
            Triangle(v31, v23, t.v3),
            Triangle(v23, v31, v12)
          )
          triangles.flatMap(subdivide(_, n-1))
        }
      }
      sphereCache(subdivisions) = sphereCache(0).flatMap(subdivide(_, subdivisions))
    }
    
    gl.glPushMatrix()
    gl.glTranslated(center.x, center.y, center.z)
    
    gl.glBegin(GL_TRIANGLES)
    for (triangle <- sphereCache(subdivisions)) {
      def sphereVertex(v: Vec3) {
        gl.glNormal3d(v.x, v.y, v.z)
        gl.glVertex3d(v.x*radius, v.y*radius, v.z*radius)
      }
      sphereVertex(triangle.v1)
      sphereVertex(triangle.v2)
      sphereVertex(triangle.v3)
    }
    gl.glEnd()
    
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
    gl.glEnable(GL_LIGHT1)
    gl.glLightfv(GL_LIGHT1, GL_AMBIENT,  Array(0.5f,0.5f,0.5f,0.2f), 0)
    gl.glLightfv(GL_LIGHT1, GL_DIFFUSE,  Array(0.9f,0.9f,0.9f,0.9f), 0)
    gl.glLightfv(GL_LIGHT1, GL_POSITION, Array(1,0.5f,1,0), 0)
  }
}
