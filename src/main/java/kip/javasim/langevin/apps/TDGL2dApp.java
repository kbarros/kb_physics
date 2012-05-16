package kip.javasim.langevin.apps;

import static scikit.util.Utilities.format;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;


public class TDGL2dApp extends Simulation {
	Grid grid = new Grid("Grid");
	TDGL2d sim;
	
	public static void main(String[] args) {
		new Control(new TDGL2dApp(), "Allen Cahn Coarsening");
	}

	public void load(Control c) {
		c.frame(grid);
		params.addm("dt", 1.0);
		params.add("L", 100.0);
		params.add("Random seed", 0);
		params.add("dx", 1.0);
		params.add("Time");
		params.add("Energy");
	}
	
	public void animate() {
		sim.readParams(params);
		
		grid.setScale(-1, 1);
		int Lp = sim.Lp;
		grid.registerData(Lp, Lp, sim.phi);

		params.set("Time", format(sim.t));
		params.set("Energy", format(sim.freeEnergy()));
	}
	
	public void clear() {
		grid.clear();
	}
	
	public void run() {
		sim = new TDGL2d(params);
		// sim.randomize();
		Job.animate();
		
		while (true) {
			for (int i = 0; i < 5; i++)
				sim.simulate();
			Job.animate();
		}
	}
}
