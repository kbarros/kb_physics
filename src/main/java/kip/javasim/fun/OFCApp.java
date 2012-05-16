package kip.javasim.fun;

import java.awt.Color;

import kip.javasim.LatticeNeighbors;
import scikit.dataset.Histogram;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;

public class OFCApp extends Simulation {
	int L;
	double dF, Fc, Fr, alpha, range;
	Plot plot = new Plot("Distribution of Earthquake Energies");
	Grid grid = new Grid("Stress Lattice");
	Histogram histogram;
	double[] stress;
	LatticeNeighbors neighbors;

	public static void main(String args[]) {
		new Control(new OFCApp(), "OFC Model");
	}

	public void load(Control c) {
		c.frame(plot, grid);
		plot.setAutoScale(true);
		params.add("Lattice size", 128);
		params.add("Critical stress (F_c)", 4.0);
		params.addm("Bin width", 100);
		params.addm("Interaction radius (R)", 4.0);
		params.addm("Loading stress (\u2206F)", 0.001);
		params.addm("Residual stress (F_R)", 1.0);
		params.addm("Dissipation (\u03B1)", 0.2);
	}

	public void animate() {
		grid.registerData(L, L, stress);
		grid.setScale(0, Fc);
		plot.registerBars("dist", new Histogram(histogram, params.iget("Bin width")), Color.RED);
		
		if (range != params.fget("Interaction radius (R)")) {
			range = params.fget("Interaction radius (R)");
			neighbors = new LatticeNeighbors(L, L, 0, range, LatticeNeighbors.Type.BORDERED);
		}
		dF = params.fget("Loading stress (\u2206F)");
		Fr = params.fget("Residual stress (F_R)");
		alpha = params.fget("Dissipation (\u03B1)");
	}
	
	public void clear() {
		plot.clear();
		grid.clear();
	}
	
	public void run() {
		L = params.iget("Lattice size");
		Fc = params.fget("Critical stress (F_c)");
		range = params.fget("Interaction radius (R)");
		neighbors = new LatticeNeighbors(L, L, 0, range, LatticeNeighbors.Type.BORDERED);
		histogram = new Histogram(1);
		histogram.setNormalizing(true);
		
		stress = new double[L*L];
		for (int i = 0; i < L*L; i++)
			stress[i] = Math.random();
		
		while (true) {
			int count = 0;
			int ndirty = 0;
			int dirty[] = new int[L*L];

			for (int i = 0; i < L*L; i++) {
				stress[i] += dF;
				if (stress[i] > Fc)
					dirty[ndirty++] = i;
			}

			while (ndirty > 0) {
				count++;
				int i = dirty[--ndirty];
				assert stress[i] > Fc;

				int[] n = neighbors.get(i);
				double unload = (stress[i] - Fr) * (1 - alpha) / n.length;
				stress[i] = Fr;

				for (int iter = 0; iter < n.length; iter++) {
					int j = n[iter];
					if (j != i) {
						if (stress[j] <= Fc && stress[j] + unload > Fc)
							dirty[ndirty++] = j;
						stress[j] += unload;
					}
				}
			}
			
			if (count > 0)
				histogram.accum(count);
			Job.animate();
		}
	}
}
