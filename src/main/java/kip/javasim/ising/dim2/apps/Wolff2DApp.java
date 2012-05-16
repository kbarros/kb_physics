package kip.javasim.ising.dim2.apps;

import java.awt.Color;

import kip.javasim.Complex;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DoubleValue;

class Ising {
	public static final double Tc = 2.0/Math.log(1.0+Math.sqrt(2.0));
	public int spin[];
	public int L;
	public int N;
	public double T;
	int[] stack;
	boolean[] cluster;
	
	
	public Ising(int _L, double _T) {
		L=_L;
		T = _T;
		N = L*L;
		spin = new int[N];
		stack = new int[N];
		cluster = new boolean[N]; 		
		for (int i = 0; i < N; i++)
			spin[i] = 1;
	}
	
	
	int neighbor(int i, int k) {
		int y = i/L;
		int x = i%L;
		int yp = (y + (k-1)%2 + L) % L;    // (k-1)%2 == {-1,  0, 1, 0}
		int xp = (x + (k-2)%2 + L) % L;    // (k-2)%2 == { 0, -1, 0, 1}
		return yp*L+xp;
	}
	
	
	public void step() {
		int nflipped = 0;
		while (nflipped < N/2) {
			double p = 1 - Math.exp(-2/T);

			stack[0] = (int)(Math.random()*N);
			cluster[stack[0]] = true;
			int stack_len = 1;

			while (stack_len-- > 0) {
				int i = stack[stack_len];
				for (int k = 0; k < 4; k++) {
					int ip = neighbor(i, k);
					if (!cluster[ip] && spin[ip] == spin[i] && Math.random() <= p) {
						cluster[ip] = true;
						stack[stack_len++] = ip;
					}
				}
				spin[i] *= -1;
				cluster[i] = false; // cannot be readded to the cluster since spin is now misaligned
				nflipped++;
			}
		}
		
		// if more than quarter of the spins have been flipped, then flip all system spins
		// for visual continuity
		int up = 0;
		for (int i = 0; i < N; i++)
			up += (spin[i]+1)/2;
		if (up < N/2) {
			for (int i = 0; i < N; i++)
				spin[i] *= -1;
		}
	}
}

public class Wolff2DApp extends Simulation {
    Grid grid = new Grid("Ising spins");
    Grid grid2 = new Grid("Conformal mapping");    
	Ising sim;
	int L, Lc;
	int original[], conformal[];
	
	public static void main(String[] args) {
		new Control(new Wolff2DApp(), "Clump Model");
	}

	public void load(Control c) {
		c.frameTogether("Spins", grid, grid2);
		params.add("L", 512);
		//params.addm("T", Ising.Tc);
		params.addm("T", new DoubleValue(Ising.Tc, Ising.Tc-0.2, Ising.Tc+0.2).withSlider());		
	}
	
	public void animate() {
		sim.T = params.fget("T");	
		ColorPalette palette = new ColorPalette();
		palette.setColor(0, new Color(80, 80, 80)); // background
		palette.setColor(-1, Color.BLACK); // spin down
		palette.setColor(1, Color.WHITE); // spin up
		palette.setColor(2, new Color(130, 200, 200)); // spin up (alternate)
		palette.setColor(4, new Color(0, 100, 0)); // branch cut
		grid.setColors(palette);
		grid2.setColors(palette);
		grid.registerData(L, L, original);
		grid2.registerData(Lc, Lc, conformal);
	}
	
	public void clear() {
		grid.clear();
		grid2.clear();
	}
	
	public void run() {
		sim = new Ising(params.iget("L"), params.fget("T"));	
    	L = sim.L;
    	Lc = L/2;
        original  = new int[L*L];
		conformal = new int[Lc*Lc];
		
        while (true) {
        	sim.step();
    		
            for (int y = 0; y < L; y++) {
            	for (int x = 0; x < L; x++) {
            		original[y*L+x] = sim.spin[y*L+x];
        			if ((x/32+y/32)%2 == 0) {
        				if (original[y*L+x] == 1)
        					original[y*L+x] = 2;
        			}
            	}
            }
            for (int yc = 0; yc < Lc; yc++) {
            	for (int xc = 0; xc < Lc; xc++) {
            		Complex c = new Complex(xc-Lc/2,yc-Lc/2);
            		c = c.sqrt().times(new Complex(1.01*L/Math.sqrt(Lc),0));
            		int x = (int)c.re() + L/2;
            		int y = (int)c.im() + L/2;
            		if (x >= 0 && x < L && y >= 0 && y < L) {
            			conformal[yc*Lc+xc+Lc/4] = original[y*L+x];
            			if (xc-Lc/2 < 0 && yc-Lc/2 > -2 && yc-Lc/2 < 2)
            				conformal[yc*Lc+xc+Lc/4] = 4; // branch cut
            		}
            		x = -(int)c.re() + L/2;
            		y = -(int)c.im() + L/2;
            		if (x >= 0 && x < L && y >= 0 && y < L) {
            			conformal[yc*Lc+xc-Lc/4] = original[y*L+x];
            			if (xc-Lc/2 < 0 && yc-Lc/2 > -2 && yc-Lc/2 < 2)
            				conformal[yc*Lc+xc-Lc/4] = 4; // branch cut
            		}
            		
            	}
           	}
            
        	Job.animate();
        }
	}
}
