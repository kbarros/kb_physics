package kip.javasim.percolation.apps;

import java.awt.Color;

import kip.javasim.Random;
import kip.javasim.ising.NewmanZiff;
import scikit.dataset.Accumulator;
import scikit.dataset.Bin;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;

public class MeanFieldPercolationApp extends Simulation {
	Plot fieldPlot = new Plot("Fields");
	int[] sites;
	NewmanZiff nz;
	Accumulator siteProfile;
	Accumulator clusterProfile;
	int L, R;
	int randomSeed;
	double p; // bond probability
	Bin maxClusterSizeBin;
	Random random;
	
	public static void main(String[] args) {
		new Control(new MeanFieldPercolationApp(), "Nucleation");
	}


	public void load(Control c) {
		c.frame(fieldPlot);
        params.add("Profile type", new ChoiceValue("Square1", "Square2", "Droplet"));
		params.add("System size", 1000);
		params.addm("Random seed", 0);
		params.addm("Probability modifier", 1.0);
		
		params.add("R");
		params.add("p*(2*R)");
		params.add("Max cluster size");
	}
	
	public void animate() {
		fieldPlot.registerLines("Site profile", siteProfile, Color.RED);
		fieldPlot.registerLines("Cluster profile", clusterProfile, Color.BLUE);
		
		params.set("Random seed", randomSeed);
		params.set("R", R);
		params.set("p*(2*R)", p*(2*R));
		params.set("Max cluster size", maxClusterSizeBin.average());
	}
    
	public void clear() {
		fieldPlot.clear();
	}
	
	public void run() {
		maxClusterSizeBin = new Bin();
		randomSeed = params.iget("Random seed");

		L = params.iget("System size");
		sites = new int[L];
		populateSites(params.sget("Profile type"));
		p *= params.fget("Probability modifier");
		double binWidth = R/50.;
		siteProfile = new Accumulator(binWidth);
		clusterProfile = new Accumulator(binWidth);
		
		while (true) {
			random = new Random(++randomSeed);
			
			bondSites();
			int maxCluster = getMaxCluster();
			maxClusterSizeBin.accum(nz.clusterSize(maxCluster) / Math.pow(L, 2./3.));
			
			for (int i = 0; i < L; i++) {
				siteProfile.accum(i, sites[i]);
				clusterProfile.accum(i, nz.clusterIndex(i) == maxCluster ? 1 : 0);
			}
			
			Job.animate();
		}
	}
	
	private void bondSites() {
		nz = new NewmanZiff(L);

		switch (0) {
		case 0:
			for (int i = 0; i < L; i++) {
				if (sites[i] != 0) {
					int nbonds = random.nextPoisson(R*p);
					for (int c = 0; c < nbonds; c++) {
						int j = (i + random.nextInt(R) + 1 + L) % L;
						if (sites[j] != 0)
							nz.addBond(i, j, 0, 0);
					}
				}
			}
			break;
			
		case 1:
			// throw bonds for mean field model G(n,p=1/n)
			for (int i = 0; i < L; i++) {
				int nbonds = random.nextPoisson((L-i-1.) / L);
				for (int c = 0; c < nbonds; c++) {
					int j = i + random.nextInt(L-i);
					nz.addBond(i, j, 0, 0);
				}
			}
			break;
			
		case 2:
			// throw bonds for mean field model, fixed edge number G(n,m=n/2)
			for (int b = 0; b < L/2; b++) {
				while (true) {
					int i = random.nextInt(L);
					int j = random.nextInt(L);
					if (sites[i] == 1 && sites[j] == 1 && i != j && !nz.isBonded(i,j)) {
						nz.addBond(i, j, 0, 0);
						break;
					}
				}
			}
			break;
		}

	}
	
	private int getMaxCluster() {
		int ret = 0;
		for (int i = 1; i < L; i++)
			if (nz.clusterSize(ret) < nz.clusterSize(i))
				ret = i;
		return nz.clusterIndex(ret);
	}
	
	private void populateSites(String profileType) {
		if (profileType.equals("Square1")) {
			R = L / 2;
			p = 0.5 / R;
			for (int i = 0; i < L; i++)
				sites[i] = 1;
		}
		else if (profileType.equals("Square2")) {
			R = L / 4;
			p = 2 * Math.atan(1./3.) / R;
			for (int i = 0; i < L; i++) {
				sites[i] = (L/4 < i && i < 3*L/4) ? 1 : 0;
			}
		}
		else if (profileType.equals("Droplet")) {
			R = L / 16;
			p = 1. / (2*R);
			for (int i = 0; i < L; i++) {
				double dx = (i-L/2.) / R;
				sites[i] = (random.nextDouble() < evalSech2(dx)) ? 1 : 0;
			}
		}
		else {
			throw new IllegalArgumentException();
		}
	}
	
	private double evalSech2(double dx) {
		double c = Math.cosh(dx);
		return 1 / (c*c);
	}
}
