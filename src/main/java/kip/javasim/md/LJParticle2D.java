package kip.javasim.md;

import kip.javasim.Vec3J;


public class LJParticle2D extends Particle {
	static final double epsilon=1;
	
	public void force(Particle that, Vec3J f) {
		_force(tag.pc.displacement(this, that), tag.radius + that.tag.radius, f);
	}
	public void force(Vec3J f) {
		_force(tag.pc.boundaryDistance(this), tag.radius, f);
	}
	public double potential(Particle that) {
		return _potential(tag.pc.displacement(this, that), tag.radius + that.tag.radius);
	}
	public double potential() {
		return _potential(tag.pc.boundaryDistance(this), tag.radius);
	}

	protected double _potential(Vec3J d, double sigma) {
		if (d == null) {
			return 0;
		}
		else {
			double r = d.norm();
			double R = tag.interactionRange;
			double a = sigma/r;
			double b = sigma/R;
			double a6 = a*a*a*a*a*a;
			double b6 = b*b*b*b*b*b; 
			double a12 = a6*a6;
			double b12 = b6*b6;
			return 4*epsilon*(a12 - 2*a6 - b12 + 2*b6);
		}
	}
	
	protected void _force(Vec3J d, double sigma, Vec3J f) {
		if (d == null) {
			f.x = f.y = f.z = 0;
		}
		else {
			double r = d.norm();
			double fmag = ljForce(r, sigma);
			f.x = fmag*d.x/r;
			f.y = fmag*d.y/r;
			f.z = fmag*d.z/r;
		}
	}
	
	protected double ljForce(double r, double sigma) {
		double a = sigma/r;
		double a6 = a*a*a*a*a*a;
		double a12 = a6*a6;
		return -(12/r)*4*epsilon*(a12 - a6);
	}
}
