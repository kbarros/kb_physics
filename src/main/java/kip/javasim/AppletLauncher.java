package kip.javasim;

import scikit.jobs.Job;

import javax.swing.JApplet;
import java.lang.reflect.*;


// needs to be reworked.
// use this class to launch a specific simulation as an applet. 
public class AppletLauncher extends JApplet {
	private static final long serialVersionUID = 1L;
	Job job;
	
	public void init() {
		String className = getParameter("Class");
		try {
			Class<?> c = Class.forName(className);
			Method m = c.getMethod("initApplet", new Class[] {Class.forName("javax.swing.JApplet")});
			job = (Job) m.invoke(null, new Object[] {this});
		} catch (ClassNotFoundException e) {
			System.err.println("AppletLauncher: Failed to find class " + className);
		} catch (NoSuchMethodException e) {
			System.err.println("AppletLauncher: method 'initApplet' not found");
		} catch (InvocationTargetException e) {
			System.err.println("AppletLauncher: method 'initApplet' not static");
		} catch (SecurityException e) {
			System.err.println("AppletLauncher: " + e);
		} catch (IllegalAccessException e) {
			System.err.println("AppletLauncher: " + e);
		} catch (IllegalArgumentException e) {
			System.err.println("AppletLauncher: " + e);
		}
	}
}

