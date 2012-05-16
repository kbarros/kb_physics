package kip.javasim.md.apps.shaker;

import static java.lang.Math.PI;
import static java.lang.Math.exp;
import static java.lang.Math.sqrt;

import java.awt.Color;

import kip.javasim.md.Particle;
import kip.javasim.md.ParticleContext;
import scikit.dataset.Function;
import scikit.dataset.Histogram;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;

public class VanHoveApp extends Simulation {
	Plot dist = new Plot("2*pi*r*G(r) versus r");
	SimulationTrajectory snapshots;
	ParticleContext pc;
	Histogram acc;
	
	public static void main(String[] args) {
		new Control(new VanHoveApp(), "Van Hove Analysis");
	}

	public void load(Control c) {
		c.frame(dist);
		params.add("Input directory", new DirectoryValue("/Users/kbarros/Desktop/data/binary/A=0.8 B=0.1 more"));
		params.add("Diffusion coef.", 0.1);
		params.add("t*", 5);
		params.add("Particle ID", 1);
	}

	public void animate() {
		final double sig2 = params.fget("Diffusion coef.")*params.fget("t*");
		Function gaussian = new Function(0, 5) {
			public double eval(double r) {
				return (2*PI*r)*(1/(sig2*2*PI))*exp(-(r*r)/(2*sig2));
			}
		};
		dist.registerLines("gaussian", gaussian, Color.RED);
		dist.registerLines("van-hove", acc, Color.BLACK);
	}

	public void clear() {
		dist.clear();
	}

	public void run() {
		snapshots = new SimulationTrajectory(params.sget("Input directory"));
		pc = snapshots.getContext();
		int id = params.iget("Particle ID");
		double tstar = params.fget("t*");
		

		acc = new Histogram(0.05);
		acc.setNormalizing(true);

		for (int i = 0; i < 150; i++) {
			double t2 = snapshots.endTime() - 1*i;
			double t1 = t2 - tstar;
			Particle[] ps1 = snapshots.get(t1);
			Particle[] ps2 = snapshots.get(t2);
			for (int j = 0; j < ps1.length; j++) {
				if (ps1[j].tag.id == id) {
					double r = sqrt(pc.displacement(ps1[j],ps2[j]).norm2());
					acc.accum(r, 2*PI*r);
				}
			}
			Job.animate();
		}
	}
}
