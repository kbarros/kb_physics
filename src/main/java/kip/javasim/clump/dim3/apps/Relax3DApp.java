package kip.javasim.clump.dim3.apps;

import static java.lang.Math.cos;
import static java.lang.Math.exp;
import static java.lang.Math.sqrt;
import static scikit.numerics.Math2.sqr;

import java.awt.Color;

import kip.javasim.clump.dim3.AbstractClump3D;
import kip.javasim.clump.dim3.FieldClump3D;
import scikit.dataset.Accumulator;
import scikit.graphics.dim2.Plot;
import scikit.graphics.dim3.Grid3D;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.util.DoubleArray;

public class Relax3DApp extends Simulation {
    Plot plot = new Plot("Structure factor");
    Plot maxplot = new Plot("Max");
    Grid3D grid = new Grid3D("test");
    Accumulator sf;
    Accumulator maxacc;
    FieldClump3D clump;

	public static void main(String[] args) {
		new Control(new Relax3DApp(), "Clump Model");
	}
	
	public void load(Control c) {
		c.frame(plot, maxplot, grid);
		params.addm("Zoom", new ChoiceValue("Yes", "No"));
		params.addm("T", 0.105);
		params.addm("dt", 1.0);
		params.add("Cutoff", 1.86418);
		params.add("R", 1000.0);
		params.add("L", 8000.0);
		params.add("dx", 250.0);
		params.add("kR bin-width", 0.1);
		params.add("Random seed", 0);
		params.add("Time");
		flags.add("Structure");
	}
	
	public void animate() {
		clump.readParams(params);
		params.set("Time", clump.time());
		
		int nc = clump.numColumns();
		grid.registerData(nc, nc, nc, clump.coarseGrained());
		
		if (flags.contains("Calc structure")) {
			sf.clear();
			clump.accumulateStructure(sf);
			plot.registerLines("Structure data", sf, Color.BLACK);
		}
		flags.clear();
		
		maxplot.registerLines("Max plot", maxacc, Color.BLACK);
	}
	
	public void clear() {
		maxplot.clear();
		plot.clear();
		grid.clear();
	}
	
	public void run() {
		maxacc = new Accumulator(0.1);
		
		clump = new FieldClump3D(params);
		clump.useNoiselessDynamics(true);
		seedField(params.fget("Cutoff"));
		
        sf = clump.newStructureAccumulator(params.fget("kR bin-width"));
        
        while (true) {
			clump.simulate();
			clump.accumulateStructure(sf);
			Job.animate();
			maxacc.accum(clump.time(), DoubleArray.max(clump.coarseGrained()));
		}
 	}
	
	void seedField(double cutoff) {
		int Lp = clump.numColumns();
		double dx = clump.dx;
		double R = clump.Rx;
		double density = AbstractClump3D.DENSITY;
		
  		for (int i = 0; i < Lp*Lp*Lp; i++) {
			double x = dx*(i%Lp - Lp/2);
			double y = dx*((i%(Lp*Lp))/Lp - Lp/2);
			double z = dx*((i/(Lp*Lp)) - Lp/2);
			double r = sqrt(x*x+y*y+z*z);
			double k = AbstractClump3D.KR_SP/R;
			
			// BCC (reciprocal lattice is FCC)
			double field = 0;
			field += cos(k * ( x + z) / sqrt(2));
			field += cos(k * (-x + z) / sqrt(2));
			field += cos(k * ( y + z) / sqrt(2));
			field += cos(k * (-y + z) / sqrt(2));
			field = (field > 0.9) ? (10) : -0.2;
			
			double mag = exp(-sqr(r/(R*cutoff)));
			clump.coarseGrained()[i] = density+mag*field;
  		}
	}
}
