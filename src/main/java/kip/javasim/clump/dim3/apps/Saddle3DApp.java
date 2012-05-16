package kip.javasim.clump.dim3.apps;

import static scikit.numerics.Math2.cube;
import static scikit.numerics.Math2.sqr;
import static scikit.util.Utilities.format;

import java.awt.Color;
import java.io.File;
import java.io.IOException;

import kip.javasim.clump.dim3.FieldClump3D;
import scikit.dataset.PointSet;
import scikit.graphics.dim2.Plot;
import scikit.graphics.dim3.Grid3D;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;
import scikit.util.Array3d;
import scikit.util.Commands;
import scikit.util.FileUtil;

public class Saddle3DApp extends Simulation {
	Grid3D grid = new Grid3D("Grid");
    Plot slice = new Plot("Slice");
	FieldClump3D clump;
	boolean adjustR;
	boolean findSaddle;
    
	double kRmin = 0.1;
	double kRmax = 14;
	double kRwidth = 0.1;

	
	public static void main(String[] args) {
		new Control(new Saddle3DApp(), "3D Clump Saddle Profile");
	}

	public void load(Control c) {
		c.frame(grid, slice);
		params.addm("Saddle", new ChoiceValue("Yes", "No"));
		params.addm("Adjust R", new ChoiceValue("No", "Yes"));
		params.addm("T", new DoubleValue(0.097, 0.0, 0.15).withSlider());
		params.addm("dt", 0.5);
		params.add("Packing", 0.05);
		params.add("Seed", new ChoiceValue("BCC", "Triangle", "Noise"));
		params.add("R", 1300.0);
		params.add("L", 16000.0);
		params.add("dx", 250.0);
		params.add("Random seed", 0);
		params.add("Time");
		params.add("F density");
		params.add("F / R^3");
		params.add("dF/dphi");
		params.add("Valid profile");
		params.add("Rx");
		params.add("Ry");
		params.add("Rz");
		flags.add("Res up");
		flags.add("Res down");
		flags.add("Dup");
	}
	
	public void animate() {
		if (flags.contains("Res up")) {
			clump.doubleResolution();
		}
		if (flags.contains("Res down")) {
			clump.halveResolution();
		}
		if (flags.contains("Dup")) {
			clump.duplicateAndEmbed();
		}
		flags.clear();
		
		findSaddle = params.sget("Saddle").equals("Yes");
		adjustR = params.sget("Adjust R").equals("Yes");
		
		clump.readParams(params);
		clump.packingFraction = params.fget("Packing");
		
		int Lp = clump.numColumns();
		grid.registerData(Lp, Lp, Lp, clump.coarseGrained());
		
		double[] sliceData = new double[Lp];
		System.arraycopy(clump.coarseGrained(), (Lp*Lp)*(Lp/2)+Lp*(Lp/2), sliceData, 0, Lp);
		slice.registerLines("Slice", new PointSet(0, 1, sliceData), Color.BLACK);
		
		double f = clump.freeEnergyDensity;
		double F_R3 = f*cube(clump.L)/(clump.Rx*clump.Ry*clump.Rz);
		params.set("dx", clump.dx);
		params.set("Rx", format(clump.Rx));
		params.set("Ry", format(clump.Ry));
		params.set("Rz", format(clump.Rz));
		params.set("Time", format(clump.time()));
		params.set("F density", format(f));
		params.set("F / R^3", format(F_R3));
		params.set("dF/dphi", format(clump.rms_dF_dphi));
		params.set("Valid profile", !clump.rescaleClipped);
	}
	
	public void clear() {
		grid.clear();
		slice.clear();
	}
	
	public void run() {
		clump = new FieldClump3D(params);
		clump.initializeFieldWithSeed(params.sget("Seed"));
		clump.useNoiselessDynamics(true);

		Job.animate();
		
		writeField();
		while (true) {
			double var1 = clump.phiVariance();
			clump.simulate();
			double var2 = clump.phiVariance();
			double scale = var1/var2;
			if (findSaddle)
				clump.scaleField(scale);
			if (adjustR) {
				clump.Rx -= clump.dt*sqr(clump.Rx)*clump.dFdensity_dRx();
				clump.Ry -= clump.dt*sqr(clump.Ry)*clump.dFdensity_dRy();
				clump.Rz -= clump.dt*sqr(clump.Rz)*clump.dFdensity_dRz();
			}
			writeField();
			Job.animate();
		}
	}
	
	double writeEndTime = 4000;
	double writeDelay = 19.9;
	double lastWrite;
	File writeDir;
	
	public void writeFiles() {
		try {
			writeDir = FileUtil.directoryDialog(null);
			System.out.println(writeDir);
			if (writeDir != null) {
				FileUtil.dumpString(writeDir+File.separator+"parameters.txt", params.toString());
				lastWrite = Double.NEGATIVE_INFINITY;
			}
		}
		catch (IOException e) {
			System.out.println(e.getMessage());
		}
	}
	
	public void loadField() {
		Array3d a = Commands.load3d();
		int len = clump.numColumns();
		System.arraycopy(a.array(), 0, clump.coarseGrained(), 0, len*len*len);
	}
	
	private void writeField() {
		if (writeDir != null && (clump.time() - lastWrite > writeDelay) && clump.time() < writeEndTime) {
			lastWrite = clump.time();
			String fname = writeDir+File.separator+"t="+format(clump.time());
			int Lp = clump.numColumns();
			Array3d a = new Array3d(Lp, Lp, Lp, clump.coarseGrained());
			a.writeFile(new File(fname));
		}
	}
}
