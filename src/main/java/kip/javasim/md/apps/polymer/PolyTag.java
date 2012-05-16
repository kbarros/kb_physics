package kip.javasim.md.apps.polymer;

import kip.javasim.md.ParticleTag;

public class PolyTag extends ParticleTag {
	PolyParticle particles[];
	
	public PolyTag(int id, PolyParticle particles[]) {
		super(id);
		this.particles = particles;
	}
	
	public PolyParticle previousParticle() {
//		int n = particles.length;
//		return particles[(id-1+n)%n];
		return (id == 0 ? null : particles[id-1]);
	}
	
	public PolyParticle nextParticle() {
		int n = particles.length;
//		return particles[(id+1)%n];
		return (id == n-1 ? null : particles[id+1]);
	}
}
