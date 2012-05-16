package kip.javasim.fun.fastestlap;

import javax.media.opengl.GL;
import javax.media.opengl.GLEventListener;
import javax.media.opengl.GLAutoDrawable;
import javax.media.opengl.glu.GLU;
import javax.media.opengl.GLCanvas;
import com.sun.opengl.util.Animator;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;
import java.awt.Cursor;
import java.awt.Frame;
import java.awt.Image;
import java.awt.Point;
import java.awt.Toolkit;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import com.sun.opengl.util.GLUT;
import java.text.NumberFormat;

class JavaRenderer implements GLEventListener, KeyListener, MouseMotionListener {
	Car car = new Car(0,0,0);
	int width, height;
	GLUT glut = new GLUT();
	static NumberFormat format = NumberFormat.getNumberInstance();;
	
	void drawString(String str) {
	     glut.glutBitmapString(GLUT.BITMAP_TIMES_ROMAN_24, str); 
	}
	
	void frameCenteredQuad(GL gl, double x, double y, double w2, double h2) {
		gl.glBegin(GL.GL_LINE_STRIP); {
			gl.glVertex2d(x-w2, y-h2);
			gl.glVertex2d(x-w2, y+h2);
			gl.glVertex2d(x+w2, y+h2);
			gl.glVertex2d(x+w2, y-h2);
			gl.glVertex2d(x-w2, y-h2);
		} gl.glEnd();
	}

	void fillCenteredQuad(GL gl, double x, double y, double w2, double h2) {
		gl.glBegin(GL.GL_QUADS); {
			gl.glVertex2d(x-w2, y-h2);
			gl.glVertex2d(x-w2, y+h2);
			gl.glVertex2d(x+w2, y+h2);
			gl.glVertex2d(x+w2, y-h2);
		} gl.glEnd();
	}

	void fillQuad(GL gl, double x0, double y0, double x1, double y1) {
		gl.glBegin(GL.GL_QUADS); {
			gl.glVertex2d(x0, y0);
			gl.glVertex2d(x0, y1);
			gl.glVertex2d(x1, y1);
			gl.glVertex2d(x1, y0);
		} gl.glEnd();
	}

	void drawReadoutBars(GL gl) {
		double BOTTOM_OFF = 0.01;
		double WIDTH = 0.8;
		double HEIGHT = 0.03;
		
		// steering
	    gl.glColor3d(0.95, 0.95, 0.95);
	    fillCenteredQuad(gl, 0, -1+HEIGHT+BOTTOM_OFF, WIDTH, HEIGHT);
	    gl.glColor3d(0.5, 0.5, 0.5);
	    frameCenteredQuad(gl, 0, -1+HEIGHT+BOTTOM_OFF, WIDTH, HEIGHT);
		double x = (WIDTH-HEIGHT) * (2*car.getSteeringPosition()-1);
		gl.glColor3d(0.1, 0.2, 0.9);
	    fillCenteredQuad(gl, x, -1+HEIGHT+BOTTOM_OFF, HEIGHT, HEIGHT);
	    
	    // acceleration/deceleration
	    gl.glColor3d(0.95, 0.95, 0.95);
	    fillCenteredQuad(gl, 0, -1+3*HEIGHT+2*BOTTOM_OFF, WIDTH/2, HEIGHT);
	    gl.glColor3d(0.5, 0.5, 0.5);
	    frameCenteredQuad(gl, 0, -1+3*HEIGHT+2*BOTTOM_OFF, WIDTH/2, HEIGHT);
	    
	    gl.glColor3d(0, 0.9, 0);
	    x = car.accel_pedal*WIDTH/2;
	    double y = -1+2*HEIGHT+2*BOTTOM_OFF;
	    fillQuad(gl, 0, y, x, y+2*HEIGHT);
	    
	    gl.glColor3d(0.9, 0.3, 0.0);
	    x = car.brake_pedal*WIDTH/2;
	    y = -1+2*HEIGHT+2*BOTTOM_OFF;
	    fillQuad(gl, 0, y, -x, y+2*HEIGHT);
	    
	}
	
	void drawOverlays(GL gl) {
		gl.glDisable(GL.GL_LIGHTING);
		gl.glDisable(GL.GL_DEPTH_TEST);
		gl.glDisable(GL.GL_BLEND);

		gl.glMatrixMode(GL.GL_PROJECTION);
		gl.glPushMatrix();
		gl.glLoadIdentity();

		gl.glOrtho(-1, 1, -1, 1, -1, 1);
		gl.glMatrixMode(GL.GL_MODELVIEW);
		gl.glPushMatrix();
		gl.glLoadIdentity();

		gl.glColor3d(0.2, 0.2, 0.2);
		gl.glRasterPos2d(-1.0, 0.9);
		drawString("" + format.format(car.km_per_hour()) + " km/h");

		drawReadoutBars(gl);

		gl.glMatrixMode(GL.GL_PROJECTION);
		gl.glPopMatrix();
		gl.glMatrixMode(GL.GL_MODELVIEW);
		gl.glPopMatrix();

		gl.glEnable(GL.GL_LIGHTING);
		gl.glEnable(GL.GL_DEPTH_TEST);
		gl.glEnable(GL.GL_BLEND);
	}
	
	public void display(GLAutoDrawable glDrawable) {
		final GL gl = glDrawable.getGL();
		gl.glClear(GL.GL_COLOR_BUFFER_BIT);
		gl.glClear(GL.GL_DEPTH_BUFFER_BIT);
		gl.glLoadIdentity();
		gl.glTranslatef(0.0f, 0.0f, -100.0f);
		
		car.simulate();
		car.draw(gl);
		drawOverlays(gl);
	}

	public void displayChanged(GLAutoDrawable gLDrawable, boolean modeChanged, boolean deviceChanged) {
	}

	public void init(GLAutoDrawable glDrawable) {
		final GL gl = glDrawable.getGL();
		gl.glShadeModel(GL.GL_FLAT);
		gl.glClearColor(1f, 1f, 1f, 0.0f);
		gl.glClearDepth(1.0f);
		gl.glDepthFunc(GL.GL_LEQUAL);
		gl.glHint(GL.GL_PERSPECTIVE_CORRECTION_HINT, GL.GL_NICEST);
		
	    gl.glLightModeli(GL.GL_LIGHT_MODEL_TWO_SIDE, GL.GL_TRUE);
	    gl.glLightfv(GL.GL_LIGHT1, GL.GL_POSITION, new float[] {1,0,1,0}, 0);
	    gl.glLightfv(GL.GL_LIGHT1, GL.GL_DIFFUSE, new float[]{1,1,1,1}, 0);
	    
		gl.glEnable(GL.GL_DEPTH_TEST);
	    gl.glEnable(GL.GL_NORMALIZE);
		gl.glEnable(GL.GL_LIGHTING);
		gl.glEnable(GL.GL_LIGHT0);
		gl.glEnable(GL.GL_LIGHT1);
		
		glDrawable.addKeyListener(this);
		glDrawable.addMouseMotionListener(this);
	}

	public void reshape(GLAutoDrawable gLDrawable, int x, int y, int width, int height) {
		this.width = width;
		this.height = height;
		
		GL gl = gLDrawable.getGL();
		gl.glMatrixMode(GL.GL_PROJECTION);
		gl.glLoadIdentity();
		(new GLU()).gluPerspective(50.0f, width/(float)height, 1.0, 1000.0);
		gl.glMatrixMode(GL.GL_MODELVIEW);
		gl.glLoadIdentity();
	}

	public void keyPressed(KeyEvent e) {
		switch (e.getKeyCode()) {
		case KeyEvent.VK_ESCAPE:
			System.exit(0);
			break;
		case KeyEvent.VK_SPACE:
			car.accelerate = true;
			break;
		case KeyEvent.VK_ALT:
			car.brake = true;
			break;
		}
	}

	public void keyReleased(KeyEvent e) {
		switch (e.getKeyCode()) {
		case KeyEvent.VK_SPACE:
			car.accelerate = false;
			break;
		case KeyEvent.VK_ALT:
			car.brake = false;
			break;
		}
	}

	public void keyTyped(KeyEvent e) {
	}

	public void mouseDragged(MouseEvent e) {
	}

	public void mouseMoved(MouseEvent e) {
		car.setSteeringPosition(e.getX()/(double)width);
	}
}


public class FastestLap implements Runnable {
	static Thread displayT = new Thread(new FastestLap());

	public static void main(String[] args) {
		displayT.start();
	}

	public void run() {
		Frame frame = new Frame("Jogl 3d Shape/Rotation");
		GLCanvas canvas = new GLCanvas();
		canvas.addGLEventListener(new JavaRenderer());
		frame.add(canvas);
		frame.setSize(1024, 768);
//		frame.setUndecorated(true);
//		int size = frame.getExtendedState();
//		size |= Frame.MAXIMIZED_BOTH;
//		frame.setExtendedState(size);
		
		frame.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				System.exit(0);
			}
		});
		frame.setVisible(true);
		canvas.requestFocus();
		
		// hide cursor
		Image cursorImage = Toolkit.getDefaultToolkit().getImage("");
		Cursor blankCursor = Toolkit.getDefaultToolkit().createCustomCursor(cursorImage, new Point( 0, 0), "" );
		canvas.setCursor( blankCursor );
		
	    Animator anim = new Animator(canvas);
	    anim.setRunAsFastAsPossible(false);
	    anim.start();
	}
}
