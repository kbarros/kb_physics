package kip.javasim.fun.fastestlap;

import javax.media.opengl.GL;
import javax.media.opengl.glu.GLU;
import javax.media.opengl.glu.GLUquadric;

import kip.javasim.Vec3J;

import com.sun.opengl.impl.GLUquadricImpl;
import static java.lang.Math.*;
import static scikit.numerics.Math2.*;

public class Car {
	final GLU glu = new GLU();
	final GLUquadric gluQuad = new GLUquadricImpl();
	
	final double body_length = 2; // meters
	final double body_width = 0.4*body_length;
	final double helmet_radius = 0.4*body_width;
	
	final double wheel_radius = 0.2 * body_length;
	final double wheel_width = 1.3 * wheel_radius;
	final double wheel_offset_w = body_width * 0.6 + wheel_width/2;
	final double wheel_offset_l = 0.4 * body_length;
	
	final double foil_distance = wheel_offset_l + 1.5*wheel_radius;
	final double foil_height = body_width * 0.6;
	final double foil_length = 2*wheel_offset_w + 0.7*wheel_width;
	final double foil_width = body_length * 0.3;
	
	final double MAX_STEER = PI/12;
	final double ENGINE_POWER = 35000; // watts, corresponds to about 47 hp
	final double MAX_ENGINE_FORCE = ENGINE_POWER / 10; // power / (36 km/hr)
	final double BRAKE_FORCE = MAX_ENGINE_FORCE;
	final double MASS = 600; // kg, about 1300 pounds
	final double INERTIA = MASS * sqr(body_length/2) / 2; // m R^2 / 2
	
	final double gravity = 9.8; // m/s^2
	final double FRICTION_S = 2;
	final double FRICTION_K = 0.6 * FRICTION_S;
	
	final double dt = 0.01;
	
	Vec3J pos;
	Vec3J vel;
	double a, va;
	
	double steer; 						// angle in radians
	boolean accelerate;
	boolean brake;
	
	final double rate_accel_pedal = 1; // seconds for accel pedal to go from 0 to 1
	final double rate_brake_pedal = 1; // seconds for decel pedal to go from 0 to 1 
	
	double accel_pedal; // between 0 and 1
	double brake_pedal; // between 0 and 1
	
	long updateMillis, lastUpdate;
	
	
	public Car(double x, double y, double angle) {
		pos = new Vec3J(x,y,0);
		vel = new Vec3J(0,0,0);
		a = angle;
		va = 0;
	}

	void drawRect(GL gl, double w, double h) {
		gl.glBegin(GL.GL_TRIANGLES); {
			gl.glVertex2d(-w/2, -h/2);
			gl.glVertex2d(-w/2, h/2);
			gl.glVertex2d(w/2, -h/2);

			gl.glVertex2d(w/2, h/2);
			gl.glVertex2d(w/2, -h/2);
			gl.glVertex2d(-w/2, h/2);
		} gl.glEnd();
	}
	
	void drawRearFoil(GL gl) {
	    gl.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE,  new float[]{0.6f, 0, 0.7f, 1}, 0);
	    gl.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_SPECULAR, new float[]{0.3f, 0.3f, 0.3f, 1}, 0);
		gl.glPushMatrix();
		gl.glTranslated(-foil_distance, 0, foil_height);
		gl.glRotated(20, 0, 1, 0);
		drawRect(gl, foil_width, foil_length);
		gl.glPopMatrix();
	}
	
	void drawBody(GL gl) {
		// draw helmet
	    gl.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE,  new float[]{0,0.2f,0.8f,1}, 0);
	    gl.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_SPECULAR, new float[]{0.3f,0.3f,0.1f,1}, 0);
		gl.glTranslated(0, 0, body_width*0.4);
		glu.gluSphere(gluQuad, helmet_radius, 10, 10);
		gl.glTranslated(0, 0, -body_width*0.4);
		
		// draw body tube
	    gl.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE,  new float[]{1,0,0,1}, 0);
	    gl.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_SPECULAR, new float[]{0.3f,0.3f,0.3f,1}, 0);
		gl.glPushMatrix();
		gl.glRotated(90, 0, 1, 0);
		gl.glTranslated(0, 0, -body_length*0.5);
		glu.gluCylinder(gluQuad, body_width*0.5, body_width*0.3, body_length*1.2, 10, 1);
		gl.glPopMatrix();
	}
	
	void drawWheel(GL gl, boolean front, double xoff, double yoff) {
		gl.glPushMatrix(); {
			gl.glTranslated(xoff, yoff, 0);
			if (front) {
				gl.glRotated(rad2deg(steer), 0, 0, 1);
			}
		    
			gl.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE,  new float[]{0.1f, 0, 0.2f, 1}, 0);
		    gl.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_SPECULAR, new float[]{0, 0, 0, 1}, 0);
			gl.glRotated(90, 1, 0, 0); 		// rotate so that wheel axel is in the z direction
			gl.glTranslated(0, 0, -wheel_width/2);
			glu.gluCylinder(gluQuad, wheel_radius, wheel_radius, wheel_width, 10, 1);

			gl.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_DIFFUSE,  new float[]{0.9f, 0.9f, 0.9f, 1}, 0);
			gl.glMaterialfv(GL.GL_FRONT_AND_BACK, GL.GL_AMBIENT, new float[]{0.6f, 0.6f, 0.6f, 1}, 0);
			glu.gluDisk(gluQuad, 0.5*wheel_radius, wheel_radius, 10, 1);
			gl.glTranslated(0, 0, wheel_width);
			glu.gluDisk(gluQuad, 0.5*wheel_radius, wheel_radius, 10, 1);
		} gl.glPopMatrix();
	}
	
	public void draw(GL gl) {
		gl.glPushMatrix();
		gl.glTranslated(pos.x, pos.y, 0);
		gl.glRotated(rad2deg(a), 0, 0, 1);
		
		drawBody(gl);
		drawRearFoil(gl);
		
		gl.glPushMatrix();
		drawWheel(gl, true,   wheel_offset_l,  wheel_offset_w);
		drawWheel(gl, true,   wheel_offset_l, -wheel_offset_w);
		drawWheel(gl, false, -wheel_offset_l,  wheel_offset_w);
		drawWheel(gl, false, -wheel_offset_l, -wheel_offset_w);
		gl.glPopMatrix();
	}
	
	/*
	// absolute offset of wheel from center of car
	public Vec3 wheelOffset(int i) {
		double w_off = (i % 2 == 0 ? 1 : -1) * wheel_offset_w;
		double l_off = (i / 2 == 0 ? 1 : -1) * wheel_offset_l;  
		return new Vec3(l_off, w_off, 0).rotateZ(a);
	}
	
	// absolute angle of direction which wheel is pointing
	public Vec3 wheelDir(int i) {
		double angle = (i < 2) ? a+steer : a;
		return new Vec3(cos(angle), sin(angle), 0);
	}

	// absolute velocity of wheel
	public Vec3 wheelVel(int i) {
		Vec3 p = wheelOffset(i);
		return vel.plus(new Vec3(0,0,va).cross(p));
	}
	*/
	
	public Vec3J wheelOffset(int i) {
		double sgn = (i == 0) ? 1 : -1;
		return new Vec3J(sgn*wheel_offset_l, 0, 0).rotateZ(a);
	}
	
	public Vec3J wheelDir(int i) {
		double angle = (i == 0) ? a+steer : a;
		return new Vec3J(1, 0, 0).rotateZ(angle);
	}
	
	public Vec3J wheelVel(int i) {
		return vel.plus(new Vec3J(0,0,va).cross(wheelOffset(i)));
	}
	
	// engine force magnitude (in direction of wheel)
	double engineForce(int i) {
		if (i == 0) // power through rear wheels only
			return 0.0;
		else {
			double speed = wheelVel(i).norm() + 1e-8;
			return min(accel_pedal*ENGINE_POWER/speed, MAX_ENGINE_FORCE); // power = force * vel
		}
	}
	
	double sign(double x) {
		return (x < 0) ? -1 : 1;
	}
	
	// brake force magnitude (in direction of wheel)
	double brakeForce(int i) {
		return - sign(wheelDir(i).dot(wheelVel(i))) * brake_pedal*BRAKE_FORCE;
	}
	
	void updateStep() {
		if (accelerate)
			accel_pedal = min(1, accel_pedal + rate_accel_pedal*dt);
		else
			accel_pedal = max(0, accel_pedal - rate_accel_pedal*dt);
		if (brake)
			brake_pedal = min(1, brake_pedal + rate_brake_pedal*dt);
		else
			brake_pedal = max(0, brake_pedal - rate_brake_pedal*dt);
		
		Vec3J forces[] = new Vec3J[] {new Vec3J(0,0,0), new Vec3J(0,0,0)};
		
		for (int i = 0; i < 2; i++) {
			Vec3J d = wheelDir(i);
			Vec3J v = wheelVel(i);
			Vec3J p = d.cross(new Vec3J(0,0,1));
			
			// net force in direction of tire given infinite friction
			double throttle = engineForce(i) + brakeForce(i);
			forces[i] = forces[i].plus(d.scale(throttle));
			
			// static friction which would be required to prevent sliding 
			Vec3J dv = v.projectOnto(p);
			forces[i] = forces[i].plus(dv.scale(- MASS / (2 * dt)));
			
			double f_s = MASS*gravity*FRICTION_S/2; // magnitude of static friction force
			double f_k = MASS*gravity*FRICTION_K/2; // magnitude of kinetic friction force
			
			System.out.println("" + i + " force " + forces[i].norm());

			// resolve between static and kinetic traction
			if (forces[i].norm() > f_s) {
				System.out.println("...slip!");
				
				// force magnitude in direction of tire
				double f_parallel = abs(throttle);
				
				// if wheels are locked then not all of the maximum break force f_k will be transferred
				// from surface to car.
				if (abs(brakeForce(i)) > abs(engineForce(i))) {
					f_parallel = min(f_parallel, abs(v.normalize().dot(d)*f_k));
				}
				else {
					f_parallel = min(f_parallel, f_k);					
					System.out.println(f_parallel);
				}
				
				// force magnitude perpendicular to tire
				double f_perp = sqrt(sqr(f_k) - sqr(f_parallel));
				
				// fix signs of forces
				f_parallel *=  sign(throttle);
				f_perp     *= -sign(v.dot(p));
				
				forces[i] = d.scale(f_parallel).plus(p.scale(f_perp));
			}
		}
		
		Vec3J force = new Vec3J(0,0,0);
		double torque = 0;
		for (int i = 0; i < 2; i++) {
			force.x += forces[i].x;
			force.y += forces[i].y;
			torque += wheelOffset(i).cross(forces[i]).z;
		}
		
		vel.x += force.x * dt / MASS;
		vel.y += force.y * dt / MASS;
		va += torque * dt / INERTIA;
		
		pos.x += vel.x * dt;
		pos.y += vel.y * dt;
		a += va*dt;
	}
	
	public void simulate() {
		long time = System.currentTimeMillis();
		lastUpdate = Math.max(lastUpdate, time - 200);
		
		while (lastUpdate < time) {
			updateStep();
			lastUpdate += Math.rint(dt * 1000);
		}
	}
	
	// mouse a double value between 0 and 1
	public void setSteeringPosition(double mouse) {
		steer = MAX_STEER * (1-2*mouse);
	}
	
	public double getSteeringPosition() {
		return (1 - steer/MAX_STEER) / 2;
	}
	
	// 5000 meters = 5 km * (1000 m / km)
	// (km / hr) = (m/s) 3600 s/hr km/(1000 m) = (m/s) 3600 / 1000
	public double km_per_hour() {
		return vel.norm() * (3600/1000);
	}
	
	double rad2deg(double rad) {
		return 180 * rad / Math.PI;
	}
}
