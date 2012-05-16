package kip.javasim.md.apps.shaker;

import static java.lang.Math.*;
import java.awt.Color;
import java.io.File;

import scikit.graphics.dim2.Scene2D;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.util.FileUtil;
import static scikit.numerics.Math2.*;
import static scikit.util.Utilities.*;

import kip.javasim.md.LJParticle2D;
import kip.javasim.md.MolecularDynamics2D;
import kip.javasim.md.ParticleContext;
import kip.javasim.md.ParticleTag;


public class LennardJonesApp extends Simulation {
	Scene2D canvas = new Scene2D("Particles");
	MolecularDynamics2D<LJParticle2D> sim;
	double lastAnimate;

	public static void main(String[] args) {
		new Control(new LennardJonesApp(), "Lennard Jones Simulation");
	}
	
	public void load(Control c) {
		c.frame(canvas);
		params.add("Output directory", new DirectoryValue(""));
		params.add("Write files", new ChoiceValue("Yes", "No"));
		params.add("Topology", new ChoiceValue("Disk", "Torus"));
		params.add("Length", 70.0);
		params.add("Area fraction A", 0.8);
		params.add("Area fraction B", 0.1);
		params.add("Radius A", 1.0);
		params.add("Radius B", 0.5);
		params.add("Epsilon", 1.0);
		params.addm("dt", 0.0025);
		params.addm("Temperature", 2);
		params.addm("Bath coupling", 0.2);
		params.add("Time");
		params.add("Reduced K.E.");
	}

	public void animate() {
		sim.setStepSize(params.fget("dt"));
		sim.setTemperature(params.fget("Temperature"), params.fget("Bath coupling"));
		params.set("Time", format(sim.time()));
		params.set("Reduced K.E.", format(sim.reducedKineticEnergy()));
		canvas.setDrawables(asList(sim.pc.boundaryDw(), sim.pc.particlesDw(sim.particles)));
	}
	
	public void clear() {
		canvas.clear();
	}

	public void run() {
		double L = params.fget("Length");
		ParticleContext pc = new ParticleContext(L, ParticleContext.typeForString(params.sget("Topology")));
		double dt = params.fget("dt");		

		ParticleTag tagA = new ParticleTag(1);
		ParticleTag tagB = new ParticleTag(2);
		tagA.pc = pc;
		tagB.pc = pc;
		tagA.radius = params.fget("Radius A");
		tagB.radius = params.fget("Radius B");
		tagA.color = new Color(0f, 0f, 1f, 0.5f);
		tagB.color = new Color(0f, 1f, 0f, 0.5f);
		double particleAreaA = PI*sqr(tagA.radius);
		double particleAreaB = PI*sqr(tagB.radius);
		double DENSITY = 1;
		tagA.mass = DENSITY*particleAreaA;
		tagB.mass = DENSITY*particleAreaB;
		double range = 3*2*max(tagA.radius,tagB.radius);
		tagA.interactionRange = range;
		tagB.interactionRange = range;
		int NA = (int) (params.fget("Area fraction A")*pc.systemArea()/particleAreaA);
		int NB = (int) (params.fget("Area fraction B")*pc.systemArea()/particleAreaB);
		
		System.out.println(NA + " " + NB);
		LJParticle2D[] particles = new LJParticle2D[NA+NB];
		for (int i = 0; i < NA; i++) {
			particles[i] = new LJParticle2D();
			particles[i].tag = tagA;
		}
		for (int i = 0; i < NB; i++) {
			particles[NA+i] = new LJParticle2D();
			particles[NA+i].tag = tagB;
		}
		
		pc.layOutParticles(particles);
		sim = new MolecularDynamics2D<LJParticle2D>(dt, pc, particles);
		lastAnimate = 0;
		
		if (params.sget("Write files").equals("No")) {
			while (true) {
				maybeAnimate();
				sim.step();
			}
		}
		else {
			File dir = FileUtil.getEmptyDirectory(params.sget("Output directory"), "output");
			FileUtil.dumpString(dir+File.separator+"parameters.txt", params.toString());
			
			dt = sim.getStepSize();
			while (sim.time()+dt/2 < 50) {
				sim.step();
				maybeAnimate();
				maybeDump(10, dir, particles);
			}
			dt = 2*dt;
			sim.setStepSize(dt);
			params.set("dt", dt);
			while (sim.time()+dt/2 < 4000) {
				sim.step();
				maybeAnimate();
				maybeDump(10, dir, particles);
			}
			while (sim.time()+dt/2 < 4500) {
				sim.step();
				maybeAnimate();
				maybeDump(1, dir, particles);
			}
			while (sim.time()+dt/2 < 4550) {
				sim.step();
				maybeAnimate();
				maybeDump(0.1, dir, particles);
			}
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
