package kip.javasim.fun;

import static java.lang.Math.PI;
import kip.javasim.Random;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.numerics.fft.FFT2D;
import scikit.numerics.fn.Function2D;
import scikit.util.Utilities;

public class CoulombPairsApp  extends Simulation {
    Grid grid = new Grid("Ising spins");
    double dt;
    CoulombSim sim;
    
	public static void main(String[] args) {
		new Control(new CoulombPairsApp(), "Ising Model");
	}

	public void load(Control c) {
		c.frame(grid);
		
		params.add("Seed", 0);
		params.add("L", 16);
		params.add("T", 0.0);
		params.addm("dt", 0.1);
		params.add("time");
	}
	
	public void animate() {
		params.set("time", Utilities.format(sim.time));
		grid.registerData(sim.L, sim.L, sim.spin);
	}
	
	public void clear() {
		grid.clear();
	}

	public void run() {
		int L = params.iget("L");
		int seed = params.iget("Seed");
		sim = new CoulombSim(L, seed);
		sim.T = params.fget("T");
		dt = params.fget("dt");
		
        while (true) {
        	sim.step(dt);            
        	Job.animate();
        }
	}	
	
}


class CoulombSim {
    int L;
    double spin[];
    double phi[];
    double interaction[];
    double time;
    double T;
	Random random = new Random();
	FFT2D fft;
	
	Function2D fourierInteraction = new Function2D() {
		public double eval(double k1, double k2) {
			if (k1 == 0 && k2 == 0) return 0;
			else return 1 / (k1*k1 + k2*k2);
		}
	};
	
    CoulombSim(int L, int seed) {
    	this.L = L;
    	spin = new double[L*L];
    	phi = new double[L*L];
		fft = new FFT2D(L, L);
		fft.setLengths(L, L);
    	interaction = buildPairInteraction();
    	random.setSeed(seed);
    	
    	for (int i = 0; i < L*L; i++)
    		spin[i] = random.nextInt(2);
    }
    
    void mcsTrial() {
    	buildPhi();
    	
    	int i = random.nextInt(L*L);
    	int j = randomNeighbor(i);
    	
    	// pair annihilation or creation
    	if (spin[i] == spin[j]) {
    	}
    	// particle move
    	else {
    		double dE = 0;
    		dE += energyCostToFlip(i);
    		flipSpin(i);
    		dE += energyCostToFlip(j);
    		flipSpin(j);
    		
    		if (dE <= 0 || random.nextDouble() < Math.exp(- dE / T)) {
    			// accept change
    		}
    		else {
    			// revert change
    			flipSpin(i);
    			flipSpin(j);
    		}
    	}
    }
    
	void step(double dt) {
		for (int i = 0; i < (int)(dt*L*L); i++) {
			mcsTrial();
		}
		time += dt;
	}
	
	double energyCostToFlip(int i) {
		return -2*(1-2*spin[i])*dressedPhi(i);
	}
	
	double dressedPhi(int i) {
		return phi[i] - spin[i]*interaction[0];
	}
	
	void flipSpin(int i) {
		spin[i] = 1 - spin[i];
		buildPhi();
	}
	
	int randomNeighbor(int i) {
		int x = i % L;
		int y = i / L;
		int r = random.nextInt(4);
		switch (r) {
		case 0: return y*L + (x+1+L)%L;
		case 1: return y*L + (x-1+L)%L;
		case 2: return ((y+1+L)%L)*L + x;
		case 3: return ((y-1+L)%L)*L + x;
		}
		return -1;
	}
	
	void buildPhi() {
		fft.convolve(spin, phi, fourierInteraction);
	}
	
	void buildPhiSlow() {
		for (int i = 0; i < L*L; i++) {
			int xi = i%L;
			int yi = i/L;
			
			phi[i] = 0;
			for (int j = 0; j < L*L; j++) {
				int xj = j%L;
				int yj = j/L;
				
				int dx = (xj-xi+L)%L;
				int dy = (yj-yi+L)%L;
				phi[i] += spin[j]*interaction[dy*L+dx];
			}
		}
	}
	
	double[] buildPairInteraction() {
		double scratch[] = new double[2*L*L];
		
		for (int x2 = -L/2; x2 < L/2; x2++) {
			for (int x1 = -L/2; x1 < L/2; x1++) {
				int i = L*((x2+L)%L) + (x1+L)%L;
				double k1 = 2*PI*x1/L;
				double k2 = 2*PI*x2/L;
				scratch[i] = fourierInteraction.eval(k1, k2);
			}
		}
		fft.transform(scratch, scratch);
		
		double ret[] = new double[L*L];
		for (int i = 0; i < L*L; i++)
			ret[i] = scratch[2*i] / (L*L);
		return ret;
	}
}
