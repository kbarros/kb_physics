package kip.javasim.md;

import static java.lang.Math.PI;
import static java.lang.Math.ceil;
import static java.lang.Math.cos;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static scikit.numerics.Math2.sqr;
import static scikit.util.Utilities.periodicOffset;

import java.awt.Color;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import kip.javasim.Random;
import kip.javasim.Vec3J;
import scikit.graphics.Drawable;
import scikit.graphics.dim2.Geom2D;
import scikit.graphics.dim2.Gfx2D;
import scikit.util.Bounds;
import scikit.util.Point;
import scikit.util.Utilities;


public class ParticleContext {
	public Random rand = new Random(0);
	public enum Type { Torus2D, Disk2D };
	public double L; // system length
	public Type type;
	
	public ParticleContext(double L, Type type) {
		this.L = L;
		this.type = type;
	}
	
	public static Type typeForString(String name) {
		if (name.equals("Disk"))
			return Type.Disk2D;
		else
			return Type.Torus2D;
	}
	
	public int dim() {
		return 2;
	}
	
	public boolean periodic() {
		switch (type) {
		case Torus2D: return true;
		case Disk2D: return false;
		}
		throw new IllegalStateException();
	}
	
	public void wrap(Point p) {
		if (periodic()) {
			p.x = (p.x+L)%L;
			p.y = (p.y+L)%L;
			p.z = (p.z+L)%L;
		}
	}
	
	public double systemArea() {
		switch (type) {
		case Disk2D:
			return PI*sqr(L/2.);
		case Torus2D:
			return L*L;
		}
		throw new IllegalStateException();
	}
	
	public Bounds getBounds() {
		return new Bounds(0, L, 0, L);
	}
	
	Vec3J tempVec = new Vec3J(0,0,0);
	
	public Vec3J boundaryDistance(Point p) {
		switch (type) {
		case Disk2D:
			double boundaryRadius = L/2.;
			double xc = p.x - L/2.;
			double yc = p.y - L/2.;
			double distanceFromCenter = sqrt(xc*xc + yc*yc);
			tempVec.x = xc*(boundaryRadius/distanceFromCenter - 1);
			tempVec.y = yc*(boundaryRadius/distanceFromCenter - 1);
			tempVec.z = 0;
			return tempVec;
		default:
			return null;
		}
	}

	public Vec3J displacement(Point p1, Point p2) {
		if (periodic()) {
			tempVec.x = periodicOffset(L, p2.x-p1.x);
			tempVec.y = periodicOffset(L, p2.y-p1.y);
			tempVec.z = 0;
		}
		else {
			tempVec.x = p2.x-p1.x;
			tempVec.y = p2.y-p1.y;
			tempVec.z = 0;
		}
		return tempVec;
	}
	
	
	public void layOutParticles(Particle[] particles) {
		int N = particles.length;
		int[] indices = Utilities.integerSequence(N);
		rand.randomizeArray(indices);
		
		switch (type) {
		case Disk2D:
			// this value of R is such that the minimum distance from a particle to the wall is half
			// the minimum interparticle distance
			double R = L / (2 + sqrt(PI/N));
			for (int i = 0; i < N; i++) {
				// these expressions for radius and angle are chosen to give a spiral trajectory
				// in which the radial distance between loops has a fixed length, and the "velocity"
				// is constant. here the particle number, i, plays the role of "time".
				// it is somewhat surprising that the angle a(t) is independent of the bounding radius R
				// and total particle number N.
				double r = R * sqrt((double)(i+1) / N);
				double a = 2 * sqrt(PI * (i+1));
				particles[indices[i]].x = L/2. + r*cos(a);
				particles[indices[i]].y = L/2. + r*sin(a);
			}
			break;
		case Torus2D:
			int rootN = (int)ceil(sqrt(N));
			double dx = L/rootN;
			for (int i = 0; i < N; i++) {
				particles[indices[i]].x = (i%rootN + 0.5 + 0.01*rand.nextGaussian()) * dx;
				particles[indices[i]].y = (i/rootN + 0.5 + 0.01*rand.nextGaussian()) * dx;
			}
			break;
		}
	}
	
	
	public static void dumpParticles(String filename, Particle[] particles) {
	  	try {
	  		DataOutputStream dos = new DataOutputStream(new FileOutputStream(filename));
	  		
	  		// write context
	  		ParticleContext pc = particles[0].tag.pc;
	  		dos.writeDouble(pc.L);
	  		dos.writeBoolean(pc.type == Type.Torus2D);
	  		
	  		// write tags
	  		Set<ParticleTag> tags = new HashSet<ParticleTag>();
	  		for (Particle p : particles)
	  			tags.add(p.tag);
	  		dos.writeInt(tags.size());
	  		for (ParticleTag tag : tags) {
	  			dos.writeInt(tag.id);
	  			dos.writeDouble(tag.mass);
	  			dos.writeDouble(tag.radius);
	  			dos.writeInt(tag.color.getRed());
	  			dos.writeInt(tag.color.getGreen());
	  			dos.writeInt(tag.color.getBlue());
	  		}
	  		
	  		// write particles
	  		dos.writeInt(particles.length);
	  		for (Particle p : particles) {
	  			dos.writeInt(p.tag.id);
	  			dos.writeDouble(p.x);
	  			dos.writeDouble(p.y);
	  			dos.writeDouble(p.z);
	  		}
	  		dos.close();
	  	}
    	catch (IOException e) {}
	}
	
	public static Particle[] readParticles(String filename) {
		Particle[] ret = null;
		try {
		    DataInputStream dis = new DataInputStream(new FileInputStream(filename));
		    
		    // read context
		    double L = dis.readDouble();
		    Type type = dis.readBoolean() ? Type.Torus2D : Type.Disk2D;
		    ParticleContext pc = new ParticleContext(L, type);
		    
		    // read tags
		    Map<Integer, ParticleTag> tags = new TreeMap<Integer, ParticleTag>();
		    int ntags = dis.readInt();
		    for (int i = 0; i < ntags; i++) {
		    	ParticleTag tag = new ParticleTag(dis.readInt());
		    	tag.pc = pc;
		    	tag.mass = dis.readDouble();
		    	tag.radius = dis.readDouble();
		    	tag.color = new Color(dis.readInt(), dis.readInt(), dis.readInt(), 128);
		    	tags.put(tag.id, tag);
		    }
		    
		    // read particles
			ret = new Particle[dis.readInt()];
			for (int i = 0; i < ret.length; i++) {
				ret[i] = new Particle();
			    ret[i].tag = tags.get(dis.readInt());
			    ret[i].x = dis.readDouble();
			    ret[i].y = dis.readDouble();
			    ret[i].z = dis.readDouble();
			}
		    dis.close();
		}
		catch (IOException e) {
			System.err.println ("Unable to read from '" + filename + "'");
		}
		return ret;
	}
	
	
	public Drawable<Gfx2D> boundaryDw() {
		switch (type) {
		case Disk2D:
			return Geom2D.circle(L/2., L/2., L/2., Color.BLACK);
		case Torus2D:
			return Geom2D.rectangle(new Bounds(0., L, 0., L), Color.BLACK);
		default:
			return null;
		}
	}
	
	public Drawable<Gfx2D> particlesDw(final Particle[] particles) {
		return new Drawable<Gfx2D>() {
			public void draw(Gfx2D g) {
				for (Particle p : particles) {
					g.setColor(p.tag.color);
					g.fillCircle(p.x, p.y, p.tag.radius);
				}
			}
			public Bounds getBounds() {
				return new Bounds(0, L, 0, L);
			}			
		};
	}
	
	public Drawable<Gfx2D> particlesLinkedDw(final Particle[] particles) {
		return new Drawable<Gfx2D>() {
			public void draw(Gfx2D g) {
				for (Particle p : particles) {
					g.setColor(p.tag.color);
					g.fillCircle(p.x, p.y, p.tag.radius);
				}
				g.setColor(Color.WHITE);
				for (int i = 0; i < particles.length-1; i++) {
					double x1 = particles[i+0].x;
					double y1 = particles[i+0].y;
					double x2 = particles[i+1].x;
					double y2 = particles[i+1].y;
					double dx = x2-x1;
					double dy = y2-y1;
					if (dx*dx + dy*dy < (L/2.)*(L/2.))
						g.drawLine(x1, y1, x2, y2);
				}
			}
			public Bounds getBounds() {
				return new Bounds(0, L, 0, L);
			}			
		};
	}
	
}
