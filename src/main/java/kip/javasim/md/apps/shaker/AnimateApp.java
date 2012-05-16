package kip.javasim.md.apps.shaker;

import static scikit.util.Utilities.asList;
import static scikit.util.Utilities.format;

import java.util.ArrayList;

import kip.javasim.md.Particle;
import kip.javasim.md.ParticleContext;
import scikit.graphics.dim2.Scene2D;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.FileValue;


public class AnimateApp extends Simulation {
	Scene2D canvas = new Scene2D("Particles");
	AbstractTrajectory snapshots;
	ParticleContext pc;
	double time;
	Particle[] particles;
	
	public static void main(String[] args) {
		new Control(new AnimateApp(), "Particle Animation").getJob().throttleAnimation(true);
	}
	
	public void load(Control c) {
		c.frame(canvas);
		params.add("Input file or directory", new FileValue("/Users/kbarros/Desktop/_c001s000100tracks.gdf"));
		params.add("Data type", new ChoiceValue("Experiment", "Simulation"));
		params.add("t start", 4000.0);
		params.add("t finish", 4500.0);
		params.add("dt", 1.0);
		params.addm("t*", 30.0);
		params.addm("r*", 0.0);
		params.add("time");
	}

	public void animate() {
		params.set("time", format(time));
		canvas.setDrawables(asList(pc.particlesDw(particles), pc.boundaryDw()));
	}
	
	public void clear() {
		canvas.clear();
	}
	
	public void run() {
		String filename = params.sget("Input file or directory");
		boolean isExperiment = params.sget("Data type").equals("Experiment");
		snapshots = isExperiment ? new ExperimentTrajectory(filename) : new SimulationTrajectory(filename);
		pc = snapshots.getContext();
		
		double ti = params.fget("t start");
		double tf = params.fget("t finish");
		double dt = params.fget("dt");
		
		for (time = ti; time < tf; time += dt) {
			double tstar = params.fget("t*");
			double rstar = params.fget("r*");
			
			if (rstar == 0) {
				particles = snapshots.get(time);
			}
			else {
				ArrayList<Particle> res = new ArrayList<Particle>();
				Particle[] ps1 = snapshots.get(time - tstar);
				Particle[] ps2 = snapshots.get(time);
				for (int i = 0; i < ps1.length; i++) {
					if (pc.displacement(ps1[i],ps2[i]).norm2() > rstar*rstar)
						res.add(ps2[i]); 
				}
				particles = res.toArray(new Particle[0]);
			}
			Job.animate();
		}
		Job.animate();
	}
}
