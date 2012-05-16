package kip.javasim.ising.dim2;

import kip.javasim.Random;

public class Ising2D {
	public int spin[];
	public int L1, L2;
	public int N;
	public double T;
	public Random random = new Random();
	public double time;
	public boolean openBoundary = false;
	public double h = 0;
	
	public Ising2D(int seed, int _L1, int _L2, double _T, boolean _openBoundary) {
		random.setSeed(seed);
		L1 = _L1;
		L2 = _L2;
		T = _T;
		openBoundary = _openBoundary;
		N = L1*L2;
		time = 0;
		
		spin = new int[N];
		randomize();
	}
	
	public void randomize() {
		for (int i = 0; i < N; i++)
			spin[i] = random.nextDouble() < 0.5 ? 1 : -1;
	}
	
	protected int neighborSum(int i) {
		int x = i % L1;
		int y = i / L1;
		int acc = 0;
		
		if (openBoundary) {
			if (y < L2-1)
				acc += spin[i+L1];
			if (y > 0)
				acc += spin[i-L1];
			if (x < L1-1)
				acc += spin[i+1];
			if (x > 0)
				acc += spin[i-1];
		}
		else {
			int yp = (y+1)%L2;
			int ym = (y-1+L2)%L2;
			int xp = (x+1)%L1;
			int xm = (x-1+L1)%L1;
			acc += spin[yp*L1+x];
			acc += spin[ym*L1+x];
			acc += spin[y*L1+xp];
			acc += spin[y*L1+xm];
		}
		
		return acc;
	}
	
	private void singleStep() {
		int i = random.nextInt(N);
		double dE = 2*spin[i]*neighborSum(i) + 2*h*spin[i];
		if (dE <= 0 || (T != 0 && random.nextDouble() < Math.exp(-dE/T))) {
			spin[i] = -spin[i];
		}
	}
	
	public void step(double mcs) {
		int n = (int) Math.max(mcs * N, 0);
		for (int i = 0; i < n; i++)
			singleStep();
		time += (double)n / N;
	}
	
	public void runUntil(double t) {
		step(t - time);
	}
}
