package kip.javasim.clump.dim2.apps;


import static java.lang.Math.max;
import static java.lang.Math.min;
import static scikit.numerics.Math2.sqr;
import static scikit.util.Utilities.format;

import java.awt.Color;

import kip.javasim.clump.dim2.FieldClump2D;
import scikit.dataset.Accumulator;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.util.DoubleArray;


public class ClumpPhase2DApp extends Simulation {
	Grid grid = new Grid("Grid");
	Plot feplot = new Plot("Free energy");
	Plot relplot = new Plot("Relaxation");
	Accumulator fe_relax;
	Accumulator fe_hex;
	Accumulator mag_hex; 
	
	FieldClump2D clump;

	public static void main(String[] args) {
		new Control(new ClumpPhase2DApp(), "Clump Model Stable Phase");
	}
	
	public void load(Control c) {
		c.frame(grid, feplot, relplot);
		params.add("dt", 0.1);
		params.add("R", 1700.0);
		params.add("L", 4000.0);
		params.add("dx", 100.0);
		params.add("T", 0.0);
//		params.add("Random seed", 0);
		params.add("Time");
		params.add("F density");
		params.add("dF/dphi");
		params.add("Rx");
		params.add("Ry");
	}

	public void animate() {
		int Lp = clump.numColumns();
		grid.registerData(Lp, Lp, clump.coarseGrained());
		
		feplot.registerPoints("Hex", fe_hex, Color.RED);
//		feplot.registerPoints("Magnitude", mag_hex, Color.BLUE);
		
		relplot.registerPoints("", fe_relax, Color.BLACK);
		
		params.set("dx", clump.dx);
		params.set("Time", format(clump.time()));
		params.set("F density", format(clump.freeEnergyDensity));
		params.set("dF/dphi", format(clump.rms_dF_dphi));
		params.set("Rx", format(clump.Rx));
		params.set("Ry", format(clump.Ry));
	}

	public void clear() {
		grid.clear();
	}

	public void run() {		
		fe_relax = new Accumulator(1);
		fe_hex = new Accumulator(0.0001);
		mag_hex = new Accumulator(0.0001);
		
		clump = new FieldClump2D(params);
		clump.initializeFieldWithHexSeed();
		clump.useNoiselessDynamics(true);
		Job.animate();
		
		clump.Rx = 1700;
		clump.Ry = 1500;
		setTemperature(0.11);
		
		simulate(0.5, 500);
		
		setTemperature(0.14);
		for (int i = 0; i < 11; i++) {
			relax();
			setTemperature(clump.T + 0.001);
		}
	}
	
	void setTemperature(double T) {
		params.set("T", format(T));
		clump.T = T;
	}
	
	public void accumFE() {
		double rmax = max(clump.Rx, clump.Ry);
		double rmin = min(clump.Rx, clump.Ry);
		if (rmax / rmin > 1.15) {
			fe_hex.accum(clump.T, clump.freeEnergyDensity);
			mag_hex.accum(clump.T, DoubleArray.max(clump.coarseGrained()));
		}
		else if (rmax / rmin < 1.02)
			System.out.println("In stable phase!");
		else
			System.out.println("Weird configuration, Rmax="+rmax+" Rmin="+rmin);
	}
	
	public void relax() {
		fe_relax.clear();
		simulate(0.5, 1000);
		accumFE();
	}
	
	void simulate(double dt, double time) {
		clump.dt = dt;
		double t1 = clump.time();
		while (clump.time() - t1 < time) {
			step();
			fe_relax.accum(clump.time(), clump.freeEnergyDensity);
		}
	}
	
	public void step() {
		clump.simulate();
		clump.Rx -= 0.1*clump.dt*sqr(clump.Rx)*clump.dFdensity_dRx();
		clump.Ry -= 0.1*clump.dt*sqr(clump.Ry)*clump.dFdensity_dRy();
		Job.animate();
	}
}
