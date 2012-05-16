package kip.javasim.clump.radial;

import static java.lang.Math.PI;
import static java.lang.Math.exp;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static scikit.numerics.Math2.j0;
import static scikit.util.Utilities.format;

import java.awt.Color;

import scikit.dataset.Function;
import scikit.dataset.PointSet;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;

public class ClumpRadialApp extends Simulation {
	Plot plot = new Plot("");
	ClumpRadial clump;
	
	public static void main(String[] args) {
		new Control(new ClumpRadialApp(), "Radial Clump");
	}
	
	// dim 2: T = 0.133
	// dim 3: T = 0.0865
	public void load(Control c) {
		c.frame(plot);
		params.addm("Saddle", new ChoiceValue("Yes", "No"));
		params.add("Dimension", new ChoiceValue("3", "2"));
		params.addm("T", 0.0865);
		params.addm("dt", 0.5);
		params.add("R", 1.);
		params.add("L", 20.);
		params.add("dx", 0.02);
		params.add("eps");
		params.add("time");
	}
	
	public Function theory() {
		final int dim = params.iget("Dimension");
		final double eps = clump.reducedTemperature();
		final double sig = 1/sqrt(2);
		final double k0 = clump.KR_SP(); 
		return new Function(0.001, 10) {
			public double eval(double x) {
				if (dim == 2)
					return 1+2*sqrt(eps)*sig*k0*j0(k0*x)*exp(-x*sqrt(eps)/sig);
				else
					return 1+(8/(sqrt(3)*PI))*sqrt(eps)*sig*k0*sin(k0*x)/(k0*x)*exp(-x*sqrt(eps)/sig);
			};	
		};
	}
	
	public void animate() {
		clump.readParams(params);
		plot.registerLines("phi", new PointSet(0, clump.dx, clump.phi), Color.BLACK);
		plot.registerLines("bar", new PointSet(0, clump.dx, clump.phibar), Color.RED);
		plot.registerLines("theory", theory(), Color.BLUE);
		
		// clump.convolveWithRange(clump.phi, clump.phibar, clump.R);
		// double[] temp = new double[clump.phi.length];
		// clump.convolveWithRangeSlow(clump.phi, temp, clump.R);
		// plot.registerLines("bar2", new PointSet(0, clump.dx, temp), Color.BLUE);
		
		params.set("eps", format(clump.reducedTemperature()));
		params.set("time", format(clump.t));
	}
	
	public void clear() {
		plot.clear();
	}
	
	public void run() {
		clump = new ClumpRadial(params);
		Job.animate();
		while (true) {
			double var1 = clump.rvariance(clump.phi);
			clump.simulate();
			double var2 = clump.rvariance(clump.phi);
			double scale = var1/var2;
			if (params.sget("Saddle").equals("Yes")) {
				clump.scaleField(scale);
			}
			Job.animate();
		}
	}
}
