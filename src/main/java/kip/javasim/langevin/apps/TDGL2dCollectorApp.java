package kip.javasim.langevin.apps;

import static scikit.util.Utilities.format;

import java.awt.Color;

import kip.javasim.ising.PercolationSite2d;

import scikit.dataset.Histogram;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;

public class TDGL2dCollectorApp extends Simulation {
	Grid grid = new Grid("Grid");
	Plot timesPlot = new Plot("Times");
	Histogram times = new Histogram(1.0);
	TDGL2d sim;
	
	public static void main(String[] args) {
		new Control(new TDGL2dCollectorApp(), "Allen Cahn Coarsening");
	}

	public void load(Control c) {
		c.frame(grid);
		c.frame(timesPlot);
		params.addm("dt", 5.0);
		params.add("L", 50.0);
		params.add("Random seed", 0);
		params.add("dx", 1.0);
		params.add("Initial cond", new ChoiceValue("Small", "Small w/ shift", "Ising"));
		params.add("Time/L^2");
		params.add("Energy");
		params.add("Stripe fraction");
	}
	
	public void animate() {		
		grid.setScale(-1, 1);
		int Lp = sim.Lp;
		grid.registerData(Lp, Lp, sim.phi);
		
		timesPlot.registerBars("Escape times", times, Color.BLUE);
		
		params.set("Time/L^2", format(sim.t));
		params.set("Energy", format(sim.freeEnergy()));
	}
	
	public void clear() {
		grid.clear();
	}
	
	public void run() {
		times = new Histogram(0.1);
		sim = new TDGL2d(params);
		
		String initialType = params.sget("Initial cond");
		
		int slabCount = 0;
		for (int iter = 1; ; iter++) {
			if ("Small".equals(initialType))
				sim.randomize();
			else if ("Small w/ shift".equals(initialType))
				sim.randomizeAndShift();
			else if ("Ising".equals(initialType))
				sim.randomizeIsing();
			else
				System.out.println("Unknown initial condition");
			
			Job.animate();
			
			
			boolean b = coarsensToStripe2();
			if (b) {
				slabCount += 1;
			}
			
			double fraction = (double)slabCount/iter;
			double plusMinus = Math.sqrt(slabCount)/iter;
			params.set("Stripe fraction", slabCount+"/"+iter+" = "+format(fraction)+" +- "+format(plusMinus));
			Job.animate();
		}
	}

	boolean coarsensToStripe1() {
		double L2 = sim.L*sim.L;
		double maxTime = 2.0*L2; // L ~ t^{1/2}
		double stripeEnergy = (2*sim.surfaceEnergyDensity*sim.L); 
					
		while (sim.t < maxTime && sim.freeEnergy() > stripeEnergy/2) {
			for (int i = 0; i < 5; i++)
				sim.simulate();
			Job.animate();
		}
		
		times.accum(sim.t / L2);
		
		return sim.freeEnergy() > stripeEnergy/2;
	}
	
	boolean coarsensToStripe2() {
		double maxTime = 400;
		while (sim.t < maxTime) {
			sim.simulate();
			Job.animate();
		}
		
		boolean openBoundary = false; // we use PBC
		PercolationSite2d nz = new PercolationSite2d(sim.Lp, sim.Lp, openBoundary);
		
		for (int i = 0; i < sim.Lp*sim.Lp; i++)
			if (sim.phi[i] > 0)
				nz.occupyAndBondSite(i);
		nz.findHomologies();
		return !nz.crossHomology() && !nz.pointHomology();
	}

}
