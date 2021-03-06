package kip.javasim.langevin.apps;

import static scikit.util.Utilities.format;
import scikit.graphics.dim3.Grid3D;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;


public class TDGL3dApp extends Simulation {
	Grid3D grid = new Grid3D("Grid");
	TDGL3d sim;
	
	public static void main(String[] args) {
		new Control(new TDGL3dApp(), "Allen Cahn Coarsening");
	}

	public void load(Control c) {
		c.frame(grid);
		params.addm("dt", 5.0);
		params.addm("r", 1.0);
		params.add("L", 15.0);
		params.add("Random seed", 41);
		params.add("dx", 0.5);
		params.add("Time");
		params.add("Energy");
	}
	
	public void animate() {
		sim.readParams(params);
		sim.r = params.fget("r");
		
		grid.setScale(-1, 1);
		int Lp = sim.Lp;
		grid.registerData(Lp, Lp, Lp, sim.phi);
		
		params.set("Time", format(sim.t));
		params.set("Energy", format(sim.freeEnergy()));
	}
	
	public void clear() {
		grid.clear();
	}
	
	public void run() {
		sim = new TDGL3d(params);
		// sim.randomize();
		Job.animate();
		
		while (true) {
			for (int i = 0; i < 1; i++)
				sim.simulate();
			Job.animate();
		}
	}
}
