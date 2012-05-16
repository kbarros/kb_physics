package kip.javasim.md.apps.shaker;

import kip.javasim.md.Particle;
import kip.javasim.md.ParticleContext;

public interface AbstractTrajectory {
	public double startTime();
	public double endTime();
	public Particle[] get(double t);	
	public ParticleContext getContext();
}
