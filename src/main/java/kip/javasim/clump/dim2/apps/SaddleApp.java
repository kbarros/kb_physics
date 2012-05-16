package kip.javasim.clump.dim2.apps;

import static scikit.util.Utilities.format;

import java.awt.Color;

import kip.javasim.clump.dim2.FieldClump2D;
import scikit.dataset.PointSet;
import scikit.graphics.GrayScale;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.util.DoubleArray;


public class SaddleApp extends Simulation {
	Grid grid = new Grid("Grid");
	Grid gmode = new Grid("Growth Mode");
	FieldClump2D clump;
	Plot plot = new Plot("Slice");
	
	public static void main(String[] args) {
		new Control(new SaddleApp(), "Clump Model Saddle Profile");
	}

	public void load(Control c) {
		c.frame(grid, plot, gmode);
		params.addm("Saddle", new ChoiceValue("Yes", "No"));
		params.addm("T", 0.14);
		params.addm("dt", 1.0);
		params.add("R", 1000.0);
		params.add("L", 30000.0);
		params.add("dx", 100.0);
		params.add("Random seed", 0);
		params.add("Time");
		params.add("F density");
		params.add("dF/dphi");
		params.add("Eigenvalue");
		params.add("del Eigenmode");
		params.add("Valid profile");
		params.add("phi max");
		flags.add("Res up");
		flags.add("Res down");
	}
	
	public void animate() {
		if (flags.contains("Res up"))
			clump.doubleResolution();
		if (flags.contains("Res down"))
			clump.halveResolution();
		flags.clear();
		clump.readParams(params);
		
		grid.setColors(new GrayScale());
		grid.setAutoScale();
//		grid.setScale(0.2, 4);
		int Lp = clump.numColumns();
		grid.registerData(Lp, Lp, clump.coarseGrained());
		
		gmode.registerData(Lp, Lp, clump.growthEigenmode);
		params.set("Eigenvalue", format(clump.growthEigenvalue));
		params.set("del Eigenmode", format(clump.rms_growthEigenmode));
		
		double[] section = new double[Lp];
		System.arraycopy(clump.growthEigenmode, Lp*(Lp/2), section, 0, Lp);
		plot.registerLines("", new PointSet(0, 1, section), Color.BLUE);
		
		params.set("dx", clump.dx);
		params.set("Time", format(clump.time()));
		params.set("F density", format(clump.freeEnergyDensity));
		params.set("dF/dphi", format(clump.rms_dF_dphi));
		params.set("phi max", format(DoubleArray.max(clump.coarseGrained())));
		params.set("Valid profile", !clump.rescaleClipped);
	}
	
	public void clear() {
		grid.clear();
		plot.clear();
	}
	
	public void run() {
		clump = new FieldClump2D(params);
		clump.useNoiselessDynamics(true);
		clump.initializeFieldWithRandomSeed();
		Job.animate();
		
		while (true) {			
			if (params.sget("Saddle").equals("Yes"))
				clump.saddleStep();
			else
				clump.simulate();
			Job.animate();
		}
	}
}
