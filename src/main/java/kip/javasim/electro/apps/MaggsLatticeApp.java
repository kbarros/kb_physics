package kip.javasim.electro.apps;

import java.awt.Color;
import java.util.Random;

import scikit.dataset.DynamicArray;
import scikit.dataset.PointSet;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import static scikit.util.DoubleArray.*;

public class MaggsLatticeApp extends Simulation {
	
	int L = 5;
	int numLinks = 3*L*L*L;
	int XDir = 0;
	int YDir = 1;
	int ZDir = 2;
	
	int x1 = 0, x2 = 1;
	double dielectric = 1.0;
	double Emag = 1.0;
	double beta = 1;
	
	Random rand = new Random(0);
	int rotatePlaquetteAccept = 0;
	int rotatePlaquetteTrials = 0;
	int shiftNetEAccept = 0;
	int shiftNetETrials = 0;
	int moveChargeAccept = 0;
	int moveChargeTrials = 0;
	
	// layout of E field is (from slowest to fastest index)
	// [z] [y] [x] [link_idx] 
	// where 'link_idx' chooses from the set of three outgoing links {XDir YDir ZDir}
	// in the "positive" direction
	double E[] = new double[L*L*L*3];
	
	DynamicArray energyData;
	Plot energyPlot = new Plot("Energy", "Time", "Histogram");
	
	
	public static void main (String[] args) {
		new Control(new MaggsLatticeApp(), "Maggs");
	}
	
	public void load(Control c) {
		c.frame(energyPlot);
		params.add("iteration");
		params.add("energy");
	}

	public void run() {
		energyData = new DynamicArray();
		initializeE();
		
		while (true) {
			for (int i = 0; i < 100; i++)
				doMCSSweep();
			energyData.append(calculateEnergy());
			Job.animate();
		}
		
	}
	
	void doMCSSweep() {
		for (int cnt = 0; cnt < numLinks; cnt++) {
			int x = rand.nextInt(L);
			int y = rand.nextInt(L);
			int z = rand.nextInt(L);
			rotatePlaquette(x, y, z);
			moveCharge(x, y, z);
		}
		shiftNetE();
	}
	
	public void animate() {
		// checkE();
		// checkAcceptanceRatios();
		
		double energyArray[] = energyData.copyArray();
		energyArray = slice(energyArray, 100, energyArray.length);
		params.set("energy", mean(energyArray)+" +- "+standardError(energyArray, 50));
		energyPlot.registerPoints("Density", new PointSet(0, 1, energyData.copyArray()), Color.BLUE);
	}
	
	public void clear() {
		energyPlot.clear();
	}
	
	
	int getLink(int x, int y, int z, int dir) {
		x = (x + L) % L;
		y = (y + L) % L;
		z = (z + L) % L;
		return (((z)*L + y)*L + x)*3 + dir;
	}
	
	void initializeE() {
		for (int i = 0; i < L*L*L*3; i++) {
			E[i] = 0;
		}
		
		// create a +1 charge at (0,0,0) and a -1 charge at (1,0,0)
		int x = 0, y = 0, z = 0;
		int i = getLink(x, y, z, XDir);
		E[i] = 1 / dielectric;
		
		// create a -1 charge at (0,1,0) and a +1 charge at (1,1,0)
//		x = 0; y = 1; z = 0;
//		i = getLink(x, y, z, XDir);
//		E[i] = -1 / dielectric;
	}
	
	int getChargeAtSite(int x, int y, int z) {
		double sum =
			E[getLink(x, y, z, XDir)] - E[getLink(x-1, y, z, XDir)] +
			E[getLink(x, y, z, YDir)] - E[getLink(x, y-1, z, YDir)] +
			E[getLink(x, y, z, ZDir)] - E[getLink(x, y, z-1, ZDir)];
		return (int)Math.round(sum*dielectric);
	}
	
	// at each site we require that $\sum_j E_{i,j} = e_i / \epsilon$ 
	void checkE() {
		for (int z = 0; z < L; z++) {
			for (int y = 0; y < L; y++) {
				for (int x = 0; x < L; x++) {
					double sum =
						E[getLink(x, y, z, XDir)] - E[getLink(x-1, y, z, XDir)] +
						E[getLink(x, y, z, YDir)] - E[getLink(x, y-1, z, YDir)] +
						E[getLink(x, y, z, ZDir)] - E[getLink(x, y, z-1, ZDir)];
					double charge = sum*dielectric;
					if (Math.abs(charge) > 1e-6) {
						System.out.println("Charge q="+charge+" at ("+x+","+y+","+z+")");
					}
				}
			}
		}
	}
	
	void checkAcceptanceRatios() {
		double r1 = rotatePlaquetteAccept/(double)rotatePlaquetteTrials;
		double r2 = shiftNetEAccept/(double)shiftNetETrials;
		double r3 = moveChargeAccept/(double)moveChargeTrials;
		System.out.println("ratios: " + r1 + " " + r2 + " " + r3);
	}
	
	double calculateEnergy() {
		double U = 0;
		for (int i = 0; i < 3*L*L*L; i++) {
			U += linkEnergy(i, 0);
		}
		return U;
	}
	
	double linkEnergy(int i, double dE) {
		return dielectric * (E[i]+dE)*(E[i]+dE) / 2.0;
	}
	
	void rotatePlaquette(int x, int y, int z) {
		
		//
		// Links on a plaquette have the structure
		//
		//      2
		//    ----->
		//   ^      ^
		// 1 |      | 3
		//   |      |
		//  (X)---->
		//      0
		//
		// where the site (x,y,z) is located at position (X)
		//
		
		// randomly select from 1 of 3 possible plaquettes
		int i0, i1, i2, i3;
		switch (rand.nextInt(3)) {
		// x-y plaquette
		case 0:
			i0 = getLink(x,y,z,XDir);
			i1 = getLink(x,y,z,YDir);
			i2 = getLink(x,y+1,z,XDir);
			i3 = getLink(x+1,y,z,YDir);
			break;
		// x-z plaquette
		case 1:
			i0 = getLink(x,y,z,XDir);
			i1 = getLink(x,y,z,ZDir);
			i2 = getLink(x,y,z+1,XDir);
			i3 = getLink(x+1,y,z,ZDir);
			break;
		// y-z plaquette
		case 2:
			i0 = getLink(x,y,z,YDir);
			i1 = getLink(x,y,z,ZDir);
			i2 = getLink(x,y,z+1,YDir);
			i3 = getLink(x,y+1,z,ZDir);
			break;
		default:
			throw new IllegalStateException();
		}
		
		// generate trial change to electric field
		double dE = Emag*(rand.nextDouble()-0.5);
		double dE0 = +dE;
		double dE1 = -dE;
		double dE2 = -dE;
		double dE3 = +dE;
		
		// calculate change in energy for the plaquette rotation
		double energy1 = linkEnergy(i0,0) + linkEnergy(i1,0) + linkEnergy(i2,0) + linkEnergy(i3,0);
		double energy2 = linkEnergy(i0,dE0) + linkEnergy(i1,dE1) + linkEnergy(i2,dE2) + linkEnergy(i3,dE3);
		double deltaEnergy = energy2 - energy1;
		
		// accept new state according to the metropolis criterion
		if (deltaEnergy < 0 || rand.nextDouble() < Math.exp(-beta*deltaEnergy)) {
			E[i0] += dE0;
			E[i1] += dE1;
			E[i2] += dE2;
			E[i3] += dE3;
			rotatePlaquetteAccept++;
		}
		rotatePlaquetteTrials++;
	}
	
	void shiftNetE() {
		double norm = Emag/Math.sqrt(numLinks);
		double dEx = norm*(rand.nextDouble()-0.5);
		double dEy = norm*(rand.nextDouble()-0.5);
		double dEz = norm*(rand.nextDouble()-0.5);
		
		double energy1 = 0;
		double energy2 = 0;
		for (int i = 0; i < L*L*L; i++) {
			// energy from x-link
			energy1 += linkEnergy(3*i+0, 0);
			energy2 += linkEnergy(3*i+0, dEx);
			// energy from y-link
			energy1 += linkEnergy(3*i+1, 0);
			energy2 += linkEnergy(3*i+1, dEy);
			// energy from z-link
			energy1 += linkEnergy(3*i+2, 0);
			energy2 += linkEnergy(3*i+2, dEz);
		}
		double deltaEnergy = energy2 - energy1;
		
		if (deltaEnergy < 0 || rand.nextDouble() < Math.exp(-beta*deltaEnergy)) {
			for (int i = 0; i < L*L*L; i++) {
				E[3*i+0] += dEx;
				E[3*i+1] += dEy;
				E[3*i+2] += dEz;
			}
			shiftNetEAccept++;
		}
		shiftNetETrials++;
	}
	
	void moveCharge(int x, int y, int z) {
		int charge1, charge2;
		int i; // link from charge1 to charge2
		
		switch (rand.nextInt(6)) {
		case 0:
			charge1 = getChargeAtSite(x-1,y,z);
			charge2 = getChargeAtSite(x, y, z);
			i = getLink(x-1,y,z,XDir);
			break;
		case 1:
			charge1 = getChargeAtSite(x, y, z);
			charge2 = getChargeAtSite(x+1,y,z);
			i = getLink(x,y,z,XDir);
			break;
		case 2:
			charge1 = getChargeAtSite(x,y-1,z);
			charge2 = getChargeAtSite(x, y, z);
			i = getLink(x,y-1,z,YDir);
			break;
		case 3:
			charge1 = getChargeAtSite(x, y, z);
			charge2 = getChargeAtSite(x,y+1,z);
			i = getLink(x,y,z,YDir);
			break;
		case 4:
			charge1 = getChargeAtSite(x,y,z-1);
			charge2 = getChargeAtSite(x, y, z);
			i = getLink(x,y,z-1,ZDir);
			break;
		case 5:
			charge1 = getChargeAtSite(x, y, z);
			charge2 = getChargeAtSite(x,y,z+1);
			i = getLink(x,y,z,ZDir);
			break;
		default: throw new IllegalStateException();
		}
		
		
		// dE is the trial change to E on the link from charge1 to charge2.
		// positive dE corresponds to _increasing_ charge1 and _decreasing_ charge2
		double dE;
		if (charge1 != 0 && charge2 == 0)
			dE = - charge1 / dielectric;
		else if (charge1 == 0 && charge2 != 0)
			dE = + charge2 / dielectric;
		else
			return;
		
		double energy1 = linkEnergy(i,0);
		double energy2 = linkEnergy(i,dE);
		double deltaEnergy = energy2 - energy1;
		
		// accept new state according to the metropolis criterion
		if (deltaEnergy < 0 || rand.nextDouble() < Math.exp(-beta*deltaEnergy)) {
			E[i] += dE;
			moveChargeAccept++;
		}
		moveChargeTrials++;
	}
}
