package kip.javasim.clump.dim1.apps;

import java.awt.Color;

import kip.javasim.clump.dim1.FieldClump1D;
import scikit.dataset.Accumulator;
import scikit.dataset.Function;
import scikit.dataset.PointSet;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;
import scikit.util.DoubleArray;


public class Clump1DApp extends Simulation {
    Plot plot = new Plot("Clump density");
    Plot sfplot = new Plot("Structure factor");
    Accumulator sf, mag;
    FieldClump1D clump;

	public static void main(String[] args) {
		new Control(new Clump1DApp(), "Clump Model");
	}
	
	public void load(Control c) {
		c.frame(plot, sfplot);
		params.addm("Noisy", new ChoiceValue("Yes", "No"));
		params.addm("T", new DoubleValue(0.09, 0, 0.3).withSlider());
		params.addm("dt", 1.0);
		params.add("R", 10000.0);
		params.add("L", 100000.0);
		params.add("dx", 100.0);
		params.add("kR bin-width", 0.1);
		params.add("Random seed", 0);
		params.add("Time");
		flags.add("Clear S.F.");
	}
	
	public void animate() {
		if (flags.contains("Clear S.F."))
			sf.clear();
		flags.clear();
		
		clump.readParams(params);
		clump.useNoiselessDynamics(!params.sget("Noisy").equals("Yes"));
		plot.registerLines("Clump data", new PointSet(0, clump.dx, clump.coarseGrained()), Color.BLACK);
		
		sfplot.registerLines("Structure data", sf, Color.BLACK);
		sfplot.registerLines("Structure theory", new Function() {
			public double eval(double kR) {
				return 1/(clump.potential(kR)/clump.T+1);
			}
		}, Color.BLUE);
	}
	
	public void clear() {
		plot.clear();
	}
	
	public void accumMagnitude() {
		mag.accum(clump.T, DoubleArray.max(clump.coarseGrained())-1);
	}
	
	public void run() {
		clump = new FieldClump1D(params);
        sf = clump.newStructureAccumulator(params.fget("kR bin-width"));
        sf.enableErrorBars(true);
        mag = new Accumulator();
        
        while (true) {
			params.set("Time", clump.time());
			clump.simulate();
			clump.accumulateStructure(sf);
			Job.animate();
		}
 	}
}

