package kip.javasim.ising.dim2.apps;

import static scikit.util.Utilities.format;
import kip.javasim.ising.dim2.PhiFourth2D;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;

public class PhiFourth2DApp extends Simulation {
	Grid grid = new Grid("Grid");
	Grid gmode = new Grid("Growth Mode");
//	Plot plot = new Plot("Slice");
	PhiFourth2D sim;
	
	public static void main(String[] args) {
		new Control(new PhiFourth2DApp(), "Phi^4 theory");
	}

	public void load(Control c) {
		c.frame(grid);
//		params.addm("Saddle", new ChoiceValue("No", "Yes"));
		params.addm("Color scale", new ChoiceValue("Absolute", "Relative"));
		params.addm("Noise", new ChoiceValue("Yes", "No"));
		params.addm("T", new DoubleValue(-0.2, -1, 1).withSlider());
		params.addm("h", new DoubleValue(0., -1, 1).withSlider());
		params.addm("anisotropy", new DoubleValue(0, -1, 1).withSlider());
		params.addm("dt", 0.3);
		params.add("R", 5.0);
		params.add("L/R", 400.0);
		params.add("dx/R", 3.125);
//		params.add("Random seed", 0);
		params.add("Time");
		params.add("F density");
//		params.add("dF/dphi");
//		params.add("Eigenvalue");
//		params.add("del Eigenmode");
//		params.add("Valid profile");
//		flags.add("Res up");
//		flags.add("Res down");
	}
	
	public void animate() {
		flags.clear();
		sim.readParams(params);
		
		if ("Absolute".equals(params.sget("Color scale"))) {
			grid.setScale(-1, 1);
		}
		else {
			grid.setAutoScale();
		}
		
		grid.setDrawRange(true);
		int Lp = sim.numColumns();
		grid.registerData(Lp, Lp, sim.phi());
		
//		gmode.registerData(Lp, Lp, clump.growthEigenmode);
//		params.set("Eigenvalue", format(clump.growthEigenvalue));
//		params.set("del Eigenmode", format(clump.rms_growthEigenmode));
		
//		double[] section = new double[Lp];
//		System.arraycopy(clump.growthEigenmode, Lp*(Lp/2), section, 0, Lp);
//		plot.registerLines("", new PointSet(0, 1, section), Color.BLUE);
		
		params.set("dx/R", sim.dx/sim.R);
		params.set("Time", format(sim.time()));
		params.set("F density", format(sim.freeEnergyDensity));
//		params.set("dF/dphi", format(sim.rms_dF_dphi));
	}
	
	public void clear() {
		grid.clear();
//		plot.clear();
	}
	
	public void run() {
		sim = new PhiFourth2D(params);
		sim.randomize();
		Job.animate();
		
		while (true) {
			for (int i = 0; i < 5; i++)
				sim.simulate();
			Job.animate();
		}
	}

}
