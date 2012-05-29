package kip.javasim.fun;

import java.awt.Color;
import java.awt.Component;

import javax.media.opengl.GL;
import javax.media.opengl.GL2;
import javax.media.opengl.GLAutoDrawable;
import javax.media.opengl.GLCapabilities;
import javax.media.opengl.GLEventListener;
import javax.media.opengl.awt.GLCanvas;
import javax.media.opengl.awt.GLJPanel;
import javax.swing.JFrame;

public class TestGfxApp {

	public static GLEventListener createListener() {
		return new GLEventListener() {
			public void display(GLAutoDrawable glDrawable) {
				GL2 gl  = glDrawable.getGL().getGL2();
				gl.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT);
				gl.glColor4fv(Color.BLUE.getComponents(null), 0);
				double x = 0, y = 0;
				double w = 100, h = 200;
				gl.glBegin(GL2.GL_QUADS);
				gl.glVertex2d(x, y);
				gl.glVertex2d(x, y+h);
				gl.glVertex2d(x+w, y+h);
				gl.glVertex2d(x+w, y);
				gl.glEnd();
			}
			public void init(GLAutoDrawable glDrawable) {
				GL gl = glDrawable.getGL();
				gl.glClearColor(1f, 1f, 1f, 0.0f);
			}
			public void reshape(GLAutoDrawable glDrawable, int x, int y, int width, int height) {
				GL gl = glDrawable.getGL();
				gl.glViewport(0, 0, width, height);
			}
			public void dispose(GLAutoDrawable glDrawable) {
			}
		};
	}
	
	public static Component createCanvas(GLEventListener listener) {
		GLCapabilities capabilities = new GLCapabilities(null);
		
		boolean useGLJPanel = true;
		if (useGLJPanel) {
			GLJPanel canvas = new GLJPanel(capabilities);
			canvas.addGLEventListener(listener);
			return canvas;
		}
		else {
			final GLCanvas canvas = new GLCanvas(capabilities);
			canvas.addGLEventListener(listener);
			return canvas;
		}
	}

	public static void main(String[] args) {
		Component component = createCanvas(createListener());
		JFrame frame = new JFrame("JOGL Test");
		frame.getContentPane().add(component);
		frame.pack();
		frame.setVisible(true);
	}
}
