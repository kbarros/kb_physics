package kip.javasim.md;

import java.awt.Color;

public class ParticleTag {
	public int id;
	public double mass, radius, charge;
	public double interactionRange;
	public Color color;
	public ParticleContext pc;
	
	public ParticleTag(int id) {
		this.id = id;
	}
}
