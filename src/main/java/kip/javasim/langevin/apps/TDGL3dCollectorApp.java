package kip.javasim.langevin.apps;

import static scikit.util.Utilities.format;

import java.awt.Color;

import scikit.dataset.Histogram;
import scikit.graphics.dim2.Plot;
import scikit.graphics.dim3.Grid3D;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;

public class TDGL3dCollectorApp extends Simulation {
	Grid3D grid = new Grid3D("Grid");
	Plot timesPlot = new Plot("Times");
	Histogram times = new Histogram(1.0);
	TDGL3d sim;
	
	public static void main(String[] args) {
		new Control(new TDGL3dCollectorApp(), "Allen Cahn Coarsening");
	}

	public void load(Control c) {
		c.frame(grid);
		c.frame(timesPlot);
		params.addm("dt", 5.0);
		params.add("L", 25.0);
		params.add("Random seed", 0);
		params.add("dx", 0.5);
		params.add("Initial cond", new ChoiceValue("Small", "Small w/ shift", "Ising"));
		params.add("Time/L^2");
		params.add("Energy");
		params.add("Slab fraction");
	}
	
	public void animate() {		
		grid.setScale(-1, 1);
		int Lp = sim.Lp;
		grid.registerData(Lp, Lp, Lp, sim.phi);
		
		timesPlot.registerBars("Escape times", times, Color.BLUE);
		
		params.set("Time/L^2", format(sim.t));
		params.set("Energy", format(sim.freeEnergy()));
	}
	
	public void clear() {
		grid.clear();
	}
	
	public void run() {
		times = new Histogram(0.1);
		sim = new TDGL3d(params);
		
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
			
			double L2 = sim.L*sim.L;
			double maxTime = 4.0*L2; // L ~ t^{1/2}
			double slabEnergy = (2*sim.surfaceEnergyDensity*L2); 
						
			while (sim.t < maxTime && sim.freeEnergy() > slabEnergy/2) {
				for (int i = 0; i < 5; i++)
					sim.simulate();
				Job.animate();
			}
			
			times.accum(sim.t / L2);
			
			if (sim.freeEnergy() > slabEnergy/2) {
				slabCount += 1;
			}
			
			double fraction = (double)slabCount/iter;
			double plusMinus = Math.sqrt(slabCount)/iter;
			params.set("Slab fraction", slabCount+"/"+iter+" = "+format(fraction)+" +- "+format(plusMinus));
			Job.animate();
		}
	}

}
