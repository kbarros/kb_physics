package kip.javasim.md.apps.shaker;

import java.awt.Color;

import kip.javasim.md.Particle;
import kip.javasim.md.ParticleContext;
import scikit.dataset.Accumulator;
import scikit.dataset.DataSet;
import scikit.dataset.DynamicArray;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;


public class AnalysisApp extends Simulation {
	Plot r2plot = new Plot("Mean squared displacement versus time");
	Plot alpha = new Plot("Non-Gaussian parameter versus time");
	SimulationTrajectory snapshots;
	ParticleContext pc;
	Accumulator dx2, dx4;
	
	public static void main(String[] args) {
		new Control(new AnalysisApp(), "Non-Gaussian Analysis");
	}
	
	public void load(Control c) {
		c.frame(r2plot, alpha);		
		params.add("Input directory", new DirectoryValue("/Users/kbarros/Desktop/data/binary/A=0.8 B=0.1 more"));
		params.addm("Log scale", new ChoiceValue("True", "False"));
		params.add("Particle ID", 1);
	}
	
	public void animate() {
		boolean ls = params.sget("Log scale").equals("True");
		r2plot.setLogScale(ls, ls);
		r2plot.registerPoints("Mean squared displacement", dx2, Color.BLUE);
		
		alpha.setLogScale(true, false);
		alpha.registerPoints("alpha", getAlpha(), Color.RED);
	}
	
	public void clear() {
		r2plot.clear();
		alpha.clear();
	}
	
	public void run() {
		snapshots = new SimulationTrajectory(params.sget("Input directory"));
		pc = snapshots.getContext();
		int id = params.iget("Particle ID");
		
		dx2 = new Accumulator(0.1);
		dx4 = new Accumulator(0.1);
		
		for (int i = 0; i < 150; i++) {
			double tf = snapshots.endTime() - 0.1*i;
			for (double time = 0.1; time < 2; time += 0.1) {
				accumMoments(tf-time, tf, id, dx2, dx4);
				Job.animate();
			}
		}
		for (int i = 0; i < 150; i++) {
			double tf = snapshots.endTime() - 1*i;
			for (double time = 1; time < 20; time += 1) {
				accumMoments(tf-time, tf, id, dx2, dx4);
				Job.animate();
			}
		}
		for (int i = 0; i < 50; i++) {
			double tf = snapshots.endTime() - 10*i;
			for (double time = 10; time < 200; time += 10) {
				accumMoments(tf-time, tf, id, dx2, dx4);
				Job.animate();
			}
		}
	}
	
	public DataSet getAlpha() {
		DynamicArray ret = new DynamicArray();		
		for (double t : dx2.keys()) {
			double d = pc.dim();
			double d2 = dx2.eval(t);
			double d4 = dx4.eval(t);
			ret.append2(t, (1/(1.+2./d)) * (d4/(d2*d2)) - 1.0);
		}
		return ret;
	}
	
	public void accumMoments(double t1, double t2, int id, Accumulator dx2, Accumulator dx4) {
		Particle[] ps1 = snapshots.get(t1);
		Particle[] ps2 = snapshots.get(t2);
		for (int i = 0; i < ps1.length; i++) {
			if (ps1[i].tag.id == id) {
				double r2 = pc.displacement(ps1[i],ps2[i]).norm2();
				dx2.accum(t2-t1, r2);
				dx4.accum(t2-t1, r2*r2);
			}
		}
	}
}



