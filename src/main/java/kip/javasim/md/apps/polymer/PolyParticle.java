package kip.javasim.md.apps.polymer;

import kip.javasim.Vec3J;
import kip.javasim.md.LJParticle2D;
import kip.javasim.md.Particle;

public class PolyParticle extends LJParticle2D {
	public void force(Particle that, Vec3J f) {
		super.force(that, f);
	}
	public void force(Vec3J f) {
		PolyTag ptag = (PolyTag)tag;
		double fx = 0, fy = 0;
		
		// interaction with boundary
		super.force(f);
		fx += f.x;
		fy += f.y;
		
		// interaction with chain neighbors
		bondForce(ptag.previousParticle(), f);
		fx += f.x;
		fy += f.y;
		bondForce(ptag.nextParticle(), f);
		fx += f.x;
		fy += f.y;
		
		// output total force
		f.x = fx;
		f.y = fy;
	}
	
	public double potential(Particle that) {
		return super.potential(that);
	}
	public double potential() {
		return super.potential();
	}
	
	public void bondForce(Particle that, Vec3J f) {
		if (that == null) {
			f.x = f.y = f.z = 0;
		}
		else {
			double sigma = tag.radius + that.tag.radius;
			Vec3J d = tag.pc.displacement(this, that);
			double r = d.norm();
			double fmag = (ljForce(r, sigma) - ljForce(2*sigma-r, sigma));
			f.x = fmag*d.x/r;
			f.y = fmag*d.y/r;
			f.z = fmag*d.z/r;
		}
	}
}
