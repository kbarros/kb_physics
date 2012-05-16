package kip.javasim.clump.dim3.apps;

import java.awt.Color;

import kip.javasim.clump.dim3.AbstractClump3D;
import kip.javasim.clump.dim3.Clump3D;
import kip.javasim.clump.dim3.FieldClump3D;
import scikit.dataset.Accumulator;
import scikit.dataset.Function;
import scikit.graphics.dim2.Plot;
import scikit.graphics.dim3.Grid3D;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;


public class Clump3DApp extends Simulation {
    Plot plot = new Plot("Structure factor");
    Grid3D grid = new Grid3D("test");
    Accumulator sf;
    AbstractClump3D clump;
    boolean fieldDynamics = false;

	public static void main(String[] args) {
		new Control(new Clump3DApp(), "Clump Model");
	}
	
	public void load(Control c) {
		c.frame(plot, grid);
		params.addm("Zoom", new ChoiceValue("Yes", "No"));
		params.addm("T", 0.09);
		params.addm("dt", 1.0);
		if (fieldDynamics) {
			params.add("R", 1000.0);
			params.add("L", 4000.0);
			params.add("dx", 250.0);
		}
		else {
			params.add("R", 4.0);
			params.add("L", 16.0);
			params.add("dx", 2.0);
		}
		params.add("kR bin-width", 0.1);
		params.add("Random seed", 0);
		params.add("Time");
		flags.add("Clear S.F.");
	}
	
	public void animate() {
		clump.readParams(params);

		if (flags.contains("Clear S.F."))
			sf.clear();
		flags.clear();
		
		int nc = clump.numColumns();
		grid.registerData(nc, nc, nc, clump.coarseGrained());
		
		plot.registerLines("Structure data", sf, Color.BLACK);
		plot.registerLines("Structure theory", new Function() {
        	public double eval(double kR) {
        		return 1/(AbstractClump3D.potential(kR)/clump.T+1);
	        }
		}, Color.BLUE);
	}
	
	public void clear() {
		plot.clear();
		grid.clear();
	}
	
	public void run() {
		clump = fieldDynamics ? new FieldClump3D(params) : new Clump3D(params);
        sf = clump.newStructureAccumulator(params.fget("kR bin-width"));
        
        while (true) {
			params.set("Time", clump.time());
			clump.simulate();
			clump.accumulateStructure(sf);
			Job.animate();
		}
 	}
}
