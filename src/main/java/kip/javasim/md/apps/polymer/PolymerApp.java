package kip.javasim.md.apps.polymer;

import static java.lang.Math.*;

import java.awt.Color;
import java.io.File;

import scikit.graphics.dim2.Scene2D;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import scikit.util.FileUtil;
import static scikit.util.Utilities.*;

import kip.javasim.md.LJParticle2D;
import kip.javasim.md.MolecularDynamics2D;
import kip.javasim.md.Particle;
import kip.javasim.md.ParticleContext;
import kip.javasim.md.ParticleTag;


public class PolymerApp extends Simulation {
	Scene2D canvas = new Scene2D("Particles");
	MolecularDynamics2D<PolyParticle> sim;
	double lastAnimate;

	public static void main(String[] args) {
		new Control(new PolymerApp(), "Polymer Simulation");
	}
	
	public void load(Control c) {
		c.frame(canvas);
		params.add("Output directory", new DirectoryValue());
		params.add("Write files", new ChoiceValue("No", "Yes"));
		params.add("Topology", new ChoiceValue("Disk", "Torus"));
		params.add("Length", 30.0);
		params.add("Particles", 120);
		params.addm("Radius", new DoubleValue(1, 0.5, 1.4).withSlider());
		params.addm("Temperature", new DoubleValue(5, 0, 10).withSlider());
		params.addm("dt", 0.005);
		params.addm("Bath coupling", 0.2);
		params.add("Time");
		params.add("Reduced K.E.");
	}

	public void animate() {
		double r = params.fget("Radius");
		for (Particle p : sim.particles)
			p.tag.radius = r;
		sim.setStepSize(params.fget("dt"));
		sim.setTemperature(params.fget("Temperature"), params.fget("Bath coupling"));
		params.set("Time", format(sim.time()));
		params.set("Reduced K.E.", format(sim.reducedKineticEnergy()));
		canvas.setDrawables(asList(sim.pc.boundaryDw(), sim.pc.particlesLinkedDw(sim.particles)));
	}
	
	public void clear() {
		canvas.clear();
	}

	public void run() {
		double L = params.fget("Length");
		double dt = params.fget("dt");		
		double radius = params.fget("Radius");
		int N = params.iget("Particles");
		PolyContext pc = new PolyContext(L, ParticleContext.typeForString(params.sget("Topology")));
		PolyParticle[] particles = new PolyParticle[N];
		
		for (int i = 0; i < N; i++) {
			ParticleTag tag = new PolyTag(i, particles);
			tag.pc = pc;
			tag.radius = radius;
			tag.color = new Color(1f-(float)i/N, 0f, (float)i/N, 0.5f);
			tag.mass = 1;
			tag.interactionRange = 6*radius;
			particles[i] = new PolyParticle();
			particles[i].tag = tag;
		}
		
		layOutParticlesSpiral(L, radius, particles);		
		sim = new MolecularDynamics2D<PolyParticle>(dt, pc, particles);
		lastAnimate = 0;
		Job.animate();
		
		if (params.sget("Write files").equals("No")) {
			while (true) {
				maybeAnimate();
				sim.step();
			}
		}
		else {
			File dir = FileUtil.getEmptyDirectory(params.sget("Output directory"), "output");
			FileUtil.dumpString(dir+File.separator+"parameters.txt", params.toString());
			while (true) {
				sim.step();
				maybeAnimate();
				maybeDump(10, dir, particles);
			}
		}
	}
	
	void layOutParticlesSpiral(double L, double radius, Particle[] particles) {
		int N = particles.length;
		double R = Math.min(L / (2 + sqrt(PI/N)), 2*radius*sqrt(N/PI));
		for (int i = 0; i < N; i++) {
			double r = R * sqrt((double)(i+1) / N);
			double a = 2 * sqrt(PI * (i+1));
			particles[i].x = L/2. + r*cos(a);
			particles[i].y = L/2. + r*sin(a);
		}		
	}
	
	void maybeAnimate() {
		long steps = round((sim.time()-lastAnimate)/sim.getStepSize()); 
		if (steps >= 10) {
			Job.animate();
			lastAnimate = sim.time();
		}
		else {
			Job.yield();
		}
	}
	
	void maybeDump(double del, File dir, LJParticle2D[] particles) {
		double dt = sim.getStepSize();
		int a = (int) ((sim.time() - dt/2)/del);
		int b = (int) ((sim.time() + dt/2)/del);
		if (b > a) {
			ParticleContext.dumpParticles(dir+File.separator+"t="+format(sim.time()), particles);			
		}
	}
}
