package kip.javasim.ising.dim2.apps;

import kip.javasim.ising.PercolationSite2d;
import kip.javasim.ising.dim2.Ising2D;
import kip.javasim.ising.dim2.IsingZero2D;
import scikit.graphics.GrayScale;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
//import scikit.util.DoubleArray;
import scikit.util.DoubleArray;
import scikit.util.Utilities;


public class Ising2DApp extends Simulation {
    Grid grid = new Grid("Ising spins");
    Grid perc = new Grid("Perc");
	Ising2D sim;
	double dt;
	double[] clusters;
	
	public static void main(String[] args) {
		new Control(new Ising2DApp(), "Ising Model");
	}

	public void load(Control c) {
		c.frame(grid, perc);
		
		params.add("Seed", 69);
		params.add("Boundary", new ChoiceValue("Periodic", "Open"));
		params.add("Root time scale", new ChoiceValue("No", "Yes"));
		params.add("L", 512);
		params.add("Ratio", 1.0);
		params.add("T", 0.0);
		params.addm("dt", 1.0);
		params.add("time");
		params.add("homology");
	}
	
	public void animate() {
		params.set("time", Utilities.format(sim.time));
		
		sim.T = params.fget("T");
		dt = params.fget("dt");
		grid.registerData(sim.L1, sim.L2, sim.spin);
		
		PercolationSite2d nz = new PercolationSite2d(sim.L1, sim.L2, sim.openBoundary);
		nz.occupyAndBondSites(sim.spin, 1);
		
		for (int i = 0; i < clusters.length; i++)
			clusters[i] = nz.clusterSize(i);
		colorSitesSpecial(clusters);
		perc.setColors(new GrayScale());
		perc.registerData(sim.L1, sim.L2, clusters);
		
		nz.findHomologies();
		params.set("homology", (nz.horizontalHomology() ? "horiz ":"")
				+ (nz.verticalHomology() ? "vert ":"")
				+ (nz.crossHomology() ? "cross ":""));
	}
	
	private void colorSitesSpecial(double a[]) {
		double m = DoubleArray.max(a);
		for (int i = 0; i < a.length; i++) {
			if (a[i] == 0)
				;
			else if (a[i] == m)
				a[i] = 2;
			else
				a[i] = 1;
		}
	}
	
//	private void arrangeSpecial() {
//		for (int i = 0; i < sim.N; i++) {
//			int x = (i % sim.L1) / (sim.L1/2);
//			int y = (i / sim.L1) / (sim.L2/2);
//			sim.spin[i] = 2*((x+y)%2) - 1;  
//		}
//	}
	
	public void clear() {
		grid.clear();
		perc.clear();
	}
	
	public void run() {
		dt = params.fget("dt");
		int seed = params.iget("Seed");
		int L2 = params.iget("L");
		int L1 = (int) (L2 * params.fget("Ratio"));
		sim = new IsingZero2D(seed, L1, L2, params.fget("T"), params.sget("Boundary").equals("Open"));
		clusters = new double[L1*L2];
		
		if (params.sget("Root time scale").equals("No")) {
			while (true) {
				Job.animate();
				sim.step(dt);
			}
		}
		else {
			double ds = dt; // reinterpret dt as ds, the scaled time step
			
			// s = 0, ds, 2ds, ...
			// t = s^2
			// dt = ((i+1)^2 - i^2) ds^2 = (2i + 1) ds^2
			for (int i = 0; ; i++) {
				Job.animate();
				sim.step((2*i+1)*ds*ds);            
			}
		}
	}
}
