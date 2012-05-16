package kip.javasim.ising.dim2.apps;

import java.awt.Color;

import kip.javasim.ising.PercolationSite2d;
import kip.javasim.ising.dim2.Ising2D;
import kip.javasim.ising.dim2.IsingZero2D;
import scikit.dataset.Accumulator;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.util.FileUtil;
import scikit.util.Utilities;


public class FreezingCollectorApp extends Simulation {
    Grid grid = new Grid("Ising spins");
    Plot homPlot = new Plot("Homology");
    Plot eventsPlot = new Plot("Topological events");
    Accumulator horz, vert, cross, topologicalEvents;
	Ising2D sim;
	double ratio;
	int count;
	
	public static void main(String[] args) {
		new Control(new FreezingCollectorApp(), "Ising Model");
	}

	public void load(Control c) {
		c.frame(grid, homPlot, eventsPlot);
		
		params.add("Write to disk", new ChoiceValue("No", "Yes"));
		params.add("Output directory", new DirectoryValue());
		params.add("Boundary", new ChoiceValue("Open", "Periodic"));
		params.add("L", 256);
		params.add("Count max", 10000);
		
		params.add("ratio");
		params.add("count");
		params.add("time");
	}
	
	public void animate() {
		params.set("ratio", ratio);
		params.set("count", count);
		params.set("time", Utilities.format(sim.time));
		
		grid.registerData(sim.L1, sim.L2, sim.spin);
		homPlot.registerPoints("Horiz", horz, Color.BLUE);
		homPlot.registerPoints("Vert", vert, Color.RED);
		homPlot.registerPoints("Cross", cross, Color.BLACK);
		eventsPlot.registerPoints("Events", topologicalEvents, Color.BLACK);
	}
	
	public void clear() {
		grid.clear();
		homPlot.clear();
		eventsPlot.clear();
	}
	
	
	public void run() {
		boolean writeToDisk = params.sget("Write to disk").equals("Yes");
		String outputDir = params.sget("Output directory");
		System.out.println(outputDir);
		String results = "";
		
		boolean openBoundary = params.sget("Boundary").equals("Open");
		int L = params.iget("L");
		
		int cntMax = params.iget("Count max");
		
//		double[] targetTimes = {50., 100., 200., 500., 2000.};
//		double[] aspectRatios = {1.};
		
		double[] targetTimes = {200.};
		double[] aspectRatios = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.};
		
		int seed = 0;
		for (double r : aspectRatios) {
			ratio = r;
			int L1 = L;
			int L2 = (int) (L/r);
		    
			horz = new Accumulator();
			vert = new Accumulator();
			cross = new Accumulator();
			horz.enableErrorBars(true);
			vert.enableErrorBars(true);
			cross.enableErrorBars(true);
			topologicalEvents = new Accumulator();
			
			for (count = 0; count < cntMax; count++) {
				// sim = new Ising2D(count, L1, L2, 0, openBoundary);
				sim = new IsingZero2D(seed++, L1, L2, 0, openBoundary);

				PercolationSite2d.Topology top1=null, top2;

				for (double t : targetTimes) {
					sim.runUntil(t);

					top2 = sampleHomology(t);
					if (top1 != null)
						topologicalEvents.accum(t, top1.equals(top2) ? 0 : 1);
					top1 = top2;

					Job.animate();
				}
			}
			
			if (writeToDisk) {
				String fname = FileUtil.getEmptyDirectory(outputDir, ""+r).toString();
				System.out.println("Writing to : " + fname);
				FileUtil.dumpColumns(fname+"/horiz.txt", horz.copyData().columns());
				FileUtil.dumpColumns(fname+"/vert.txt", vert.copyData().columns());
				FileUtil.dumpColumns(fname+"/cross.txt", cross.copyData().columns());
				results +=
					r+" "+
					horz.eval(horz.maxKey())+" "+
					horz.evalError(horz.maxKey())+" "+
					vert.eval(vert.maxKey())+" "+
					vert.evalError(vert.maxKey())+" "+
					cross.eval(cross.maxKey())+" "+
					cross.evalError(cross.maxKey())+"\n";
			}
		}
		
		if (writeToDisk) {
			String fname = outputDir + "/collected.txt";
			System.out.println("Writing to : " + fname);
			FileUtil.dumpString(fname, results);
		}
	}
	
	private PercolationSite2d.Topology sampleHomology(double targetTime) {
		PercolationSite2d nz = new PercolationSite2d(sim.L1, sim.L2, sim.openBoundary);
		nz.occupyAndBondSites(sim.spin, 1);
		nz.findHomologies();
		horz.accum(targetTime, nz.horizontalHomology() ? 1 : 0);
		vert.accum(targetTime, nz.verticalHomology() ? 1 : 0);
		cross.accum(targetTime, (nz.crossHomology() || nz.pointHomology()) ? 1 : 0);
		return nz.getTopology();
	}
}
