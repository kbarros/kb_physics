package kip.javasim.glass.apps;


import kip.javasim.Random;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DoubleValue;
import scikit.jobs.params.Parameters;
import scikit.numerics.Math2;



class SpinGlass {
	public int spin[];
	public int L;
	public int N;
	public double T;
	public int range;
	
	public double D0, h0;
	public double Dv = 0, D[];
	public double hv = 0, h[];
	Random random = new Random();
	
	
	public SpinGlass(Parameters params) {
		L = params.iget("L");
		T = params.fget("T");
		N = L*L;

		random.setSeed(params.iget("Random seed"));
		
		spin = new int[N];
		D = new double[N];
		h = new double[N];
		
		Dv = params.fget("D variance");
		hv = params.fget("h variance");
		for (int i = 0; i < N; i++) {
			spin[i] = random.nextDouble() < 0.5 ? 1 : -1;
			D[i] = Dv*random.nextGaussian();
			h[i] = hv*random.nextGaussian();
		}
	}
	
	public void readParams(Parameters params) {
		T = params.fget("T");	
		range = params.iget("Range");
		D0 = params.fget("D mean");
		h0 = params.fget("h mean");
	}
	
	int neighbor(int i, int k) {
		int y = i/L;
		int x = i%L;
		int yp = (y + (k-1)%2 + L) % L;    // (k-1)%2 == {-1,  0, 1, 0}
		int xp = (x + (k-2)%2 + L) % L;    // (k-2)%2 == { 0, -1, 0, 1}
		return yp*L+xp;
	}
	
	
	// nearest neighbor ferromagnet interaction
	double srKernel(int i) {
		double ret = 0;
		for (int k = 0; k < 4; k++) {
			int j = neighbor(i, k);
			ret += spin[j];
		}
		return 10*ret;
	}
	
	
	// long range anisotropic interaction
	double lrKernel(int i) {
		int y = i/L;
		int x = i%L;
		
		double ret = 0;
		
		for (int dx = -range; dx <= +range; dx++) {
			for (int dy = -range; dy <= +range; dy++) {
				double r2 = dx*dx + dy*dy;
				if (0 < r2 && r2 < range*range) {
					double cosTheta2 = (dy*dy) / r2; // cos(theta)^2
					double cos4Theta = 8 * cosTheta2 * (cosTheta2 - 1); // cos(4 theta)
					
					int yp = (y + dy + L) % L;
					int xp = (x + dx + L) % L;
					int j = yp*L + xp;
					
					double eta = 0.5;
					double angular = (1 - eta*cos4Theta) / Math2.sqr(eta - cos4Theta);
					ret -= angular * spin[j] / r2;
				}
			}
		}
		
		return ret;
	}
	
	public void step() {
		int i = random.nextInt(N);
		int s = spin[i];
		int sp = random.nextInt(3) - 1; // next spin, randomly selected spin {-1,0,1}
		
		double kernel = (srKernel(i) + lrKernel(i));
		double E0 = (D0+D[i])*s*s   - kernel*s  - (h0+h[i])*s;
		double E1 = (D0+D[i])*sp*sp - kernel*sp - (h0+h[i])*sp;
		double dE = E1-E0;
		if (dE <= 0 || (T != 0 && random.nextDouble() < Math.exp(-dE/T))) {
			spin[i] = sp;
		}
	}
}


public class SpinGlassApp extends Simulation {
    Grid grid = new Grid("Ising spins");
	SpinGlass sim;
	int L;
	
	public static void main(String[] args) {
		new Control(new SpinGlassApp(), "Strain glass");
	}

	public void load(Control c) {
		grid.setScale(-1, +1);
		c.frame(grid);
		params.add("L", 32);
		params.add("Random seed", 0);
		params.add("D variance", 0.0);
		params.add("h variance", 0.0);
		params.addm("D mean", new DoubleValue(0, -20, 20).withSlider());
		params.addm("h mean", new DoubleValue(0, -10, 10).withSlider());
		params.addm("T", new DoubleValue(1, 0, 20).withSlider());
		params.addm("Range", 5);
	}
	
	public void animate() {
		sim.readParams(params);
		grid.registerData(L, L, sim.spin);
	}
	
	public void clear() {
		grid.clear();
	}
	
	public void run() {
		sim = new SpinGlass(params);	
    	L = sim.L;
		
        while (true) {
        	sim.step();           
        	Job.animate();
        }
	}
}
