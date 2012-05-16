package kip.javasim.fun;

import static java.lang.Math.exp;
import static java.lang.Math.min;
import static java.lang.Math.random;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;


public class SpinodalApp extends Simulation {
    Grid grid = new Grid("Grid");
    int L;
    double[] data;
    double beta, J;
    
	public static void main(String[] args) {
		new Control(new SpinodalApp(), "Spinodal Model");
	}
    
	public void load(Control c) {
		c.frame(grid);
        params.add("Size", 128);
        params.addm("Temperature", 2.0);
        params.addm("Interaction", 0.5);
    }
    
    public void animate() {
        beta = 1 / params.fget("Temperature");
        J = params.fget("Interaction");
        grid.setScale(0, 16);
        grid.registerData(L, L, data);
    }
    
    public void clear() {
    	grid.clear();
    }
    
    private double sumNeighbors(int i) {
        int N = L*L;
        int y = i/L;
        int up   = (i+L)%N;
        int down = (i-L+N)%N;
        int left = (i-1+L)%L + y*L;
        int rght = (i+1)%L + y*L;
        return data[up] + data[down] + data[left] + data[rght];
    }
    
    private double deltaEnergy(int i, int dn) {
        return dn*(dn + 2*data[i] - J*sumNeighbors(i));
    }
    
    private void monteCarloTrial() {
        int i = (int) (L*L * random());
        int dn = random() < 0.5 ? -1 : 1;
        
        if (data[i] + dn >= 0) {
            double dE = deltaEnergy(i, dn);
            if (random() < min(exp(-beta*dE), 1)) {
                data[i] += dn;
            }
        }
    }
    
    public void run() {
        L = params.iget("Size");
        data = new double[L*L];
        
        while (true) {
            for (int i = 0; i < L*L; i++)
                data[i] = 0;
            Job.animate();
            
            while (true) {
                for (int i = 0; i < L*L; i++) {
                    monteCarloTrial();
                }
                Job.animate();
            }
        }
    }
}
