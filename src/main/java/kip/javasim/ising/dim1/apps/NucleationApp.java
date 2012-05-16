package kip.javasim.ising.dim1.apps;


import static java.lang.Math.*;
import static scikit.util.Utilities.*;

import java.awt.Color;

import kip.javasim.ising.dim1.Dynamics1D;
import kip.javasim.ising.dim1.PhiFourth;
import scikit.graphics.dim2.Geom2D;
import scikit.graphics.dim2.Plot;
import scikit.dataset.Histogram;
import scikit.dataset.Accumulator;
import scikit.dataset.Function;
import scikit.dataset.PointSet;
import scikit.jobs.*;
import scikit.jobs.params.ChoiceValue;

public class NucleationApp extends Simulation {
	Plot fieldPlot = new Plot("Fields");
	Plot profilePlot = new Plot("Average Droplet Profile");
	Plot nucTimes = new Plot("Nucleation Times"); // bw 0.1
    Histogram nucTimesAcc;
    Accumulator droplet;
    Function saddleProfile;
    
	boolean phifour = true;
	Dynamics1D sim;
	
	// range of time for which to collect nucleating droplets
	double EARLY_END = 15;
    double LATE_BEGIN = 30;
	
	public static void main(String[] args) {
		new Control(new NucleationApp(), "Nucleation");
	}
	
	public void load(Control c) {
		c.frame(fieldPlot, profilePlot, nucTimes);
		params.add("Memory time", 20.0);
        params.add("Profile type", new ChoiceValue("None", "Early", "Late"));
		params.add("Data path", "");
        params.add("Profile count", 0);
		params.addm("Max count", 50000);
		params.addm("Random seed", 0);
		params.addm("Bin width", 0.5);
		if (phifour) {
			double eps = -1;
			params.add("N/R", 300.0);
			params.add("dx/R", 1.0);
			params.addm("R", 2000);
			params.addm("dt", 0.1);
			params.addm("h", sqrt(-8*eps/27) - 0.005);
			params.addm("\u03b5", eps);
		}
	}
	
	public void animate() {
		sim.setParameters(params);
		
		double bw = params.fget("Bin width");
		nucTimes.registerBars("Nucleation times", new Histogram(nucTimesAcc, bw), Color.RED);
        
		double dx_R = (double)sim.dx/sim.R;
		double N_R = (double)sim.N/sim.R;
		fieldPlot.registerLines("Field", new PointSet(0, dx_R, sim.copyField()), Color.BLACK);
		fieldPlot.setDrawables(asList(Geom2D.line(0., 0., N_R, 0., Color.BLUE)));
        
		profilePlot.registerLines("Data", droplet, Color.RED);
		profilePlot.registerLines("Theory", saddleProfile, Color.BLUE);
	}
    
	public void clear() {
		nucTimes.clear();
		fieldPlot.clear();
		profilePlot.clear();
	}
	
    void accumulateDroplet(double[] field, double pos) {
        int j = (int)round(pos/sim.dx);
        int c = sim.N/sim.dx;
        for (int i = 0; i < field.length; i++) {
			double x = (double)(i-field.length/2)*sim.dx/sim.R;
            droplet.accum(x, field[(i+j+c/2)%c]);
        }
    }
    
	void equilibrate() {
		sim.h = -sim.h;
		sim.runUntil(10);
		sim.h = -sim.h;
		sim.resetTime();
	}
	
	void simulateUntilNucleation() {
        double lowBound = 0;
        double highBound = Double.MAX_VALUE;
        boolean findProfile = false;
        if (params.sget("Profile type").equals("Early")) {
            highBound = EARLY_END;
            findProfile = true;
        }
        else if (params.sget("Profile type").equals("Late")) {
            lowBound = LATE_BEGIN;
            findProfile = true;
        }
        
		while (!sim.nucleated() && sim.time() < highBound) {
			sim.step();
			Job.animate();
		}
        
        if (lowBound < sim.time() && sim.time() < highBound) {
            if (findProfile) {
                // average difference between crude nucleation time, and real nucleation time
                double overshootEstimate = params.fget("Memory time")/2;
                double[] drop = sim.nucleationTimeAndLocation(overshootEstimate);
                accumulateDroplet(sim.simulationAtTime(drop[0]).copyField(), drop[1]);
            }
            nucTimesAcc.accum(sim.time());
            params.set("Profile count", params.iget("Profile count")+1);            
        }
	}
	
	
	public void run() {
		if (phifour) {
			sim = new PhiFourth(params);
		}
		sim.initialize(params);
		saddleProfile = sim.clone().saddleProfile();
		
		nucTimesAcc = new Histogram(params.fget("Bin width"));
        nucTimesAcc.setNormalizing(true);
        droplet = new Accumulator((double)sim.dx/sim.R);
        
		while (params.iget("Profile count") < params.iget("Max count")) {
            sim.initialize(params);
            equilibrate();
            simulateUntilNucleation();
            params.set("Random seed", params.iget("Random seed")+1);
            Job.animate();
		}
	}
}
