package kip.javasim.ising.dim2.apps;

import static scikit.util.Utilities.format;
import kip.javasim.Random;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DoubleValue;
import scikit.jobs.params.Parameters;
import scikit.numerics.fft.FFT2D;
import scikit.numerics.fn.Function2D;

class BinaryPFC2d {
	Random random = new Random();
	double L, dx, dt, t;
	int Lp;
	
	double r; // reduced temperature
	double a0, a1; // lattice constants
	double rho0, rho1; // net densities
	double noise;
	
	double[] n0, n1; // densities
	double[] scratch1, scratch2;
	FFT2D fft;
	
	double[] potential00;
	double[] potential11;
	double[] potential01;
	
	
	public BinaryPFC2d(Parameters params) {
		random.setSeed(params.iget("Random seed", 0));
		
		L = params.fget("L");
		dx = params.fget("dx");
		Lp = (int) (L / dx);
		dx = L / Lp;
		params.set("dx", dx);
		
		n0 = new double[Lp*Lp];
		n1 = new double[Lp*Lp];
		scratch1 = new double[Lp*Lp];
		scratch2 = new double[Lp*Lp];
		
		potential00 = new double[Lp*Lp];
		potential11 = new double[Lp*Lp];
		potential01 = new double[Lp*Lp];
		
		fft = new FFT2D(Lp, Lp);
		fft.setLengths(L, L);
		
		rho0 = params.fget("rho0");
		rho1 = params.fget("rho1");
		readParams(params);
		
		double[] pot = new double[Lp];
		for (int i = 0; i < Lp; i++) pot[i] = potential00[i];
//		Commands.plot(pot);
		
		
		t = 0;
		for (int i = 0; i < Lp*Lp; i++) {
			if (i%Lp < Lp/2) {
				n0[i] = rho0 + noise();
				n1[i] = 0.1;
			}
			else {
				n0[i] = 0.1;
				n1[i] = rho1 + noise();
			}
		}
	}
	
		
	public void readParams(Parameters params) {
		r = params.fget("r");
		a0 = params.fget("a0");
		a1 = params.fget("a1");
		
		dt = params.fget("dt");
		noise = params.fget("Noise");
		
		buildPotentials();
	}
	
	public void buildPotentials() {
		fft.buildFourierArrayFromReal(potential00, new Function2D() {
			public double eval(double x, double y) {
				return (x*x+y*y < a0*a0) ? 1 : 0;
			}
		});
		
		fft.buildFourierArrayFromReal(potential11, new Function2D() {
			public double eval(double x, double y) {
				return (x*x+y*y < a1*a1) ? 1 : 0;
			}
		});
		
		fft.buildFourierArrayFromReal(potential01, new Function2D() {
			public double eval(double x, double y) {
				return (Math.sqrt(x*x+y*y) < (a0+a1)/2) ? 5 : 0;
			}
		});

		// regulator
		fft.buildFourierArray(scratch1, new Function2D() {
			public double eval(double kx, double ky) {
				return (kx*kx+ky*ky > 10*10) ? 1 : 0;
			}
		});		
		for (int i = 0; i < Lp*Lp; i++) {
			potential00[i] += scratch1[i];
			potential11[i] += scratch1[i];
			potential01[i] += scratch1[i];
		}
	}
	
	public void laplacian(double[] src, double[] dst, final double cutoff) {
		fft.convolve(src, dst, new Function2D() {
			public double eval(double kx, double ky) {
				double k2 = kx*kx + ky*ky;
				return -Math.min(k2, cutoff);
			}
		});
	}
	
	public void simulate() {
		double cpl = 0;
		double trim = 1;
		
		//
		// Update n0
		//
		
		fft.convolve(n0, scratch1, potential00);
		fft.convolve(n1, scratch2, potential01);
		for (int i = 0; i < Lp*Lp; i++) {
			scratch1[i] = r*n0[i] + scratch1[i] + cpl*scratch2[i] + n0[i]*n0[i]*n0[i] - trim*1/n0[i]; // df/dn0 = r n0 + v00 n0 + v01 n1 + n^3  
		}
		laplacian(scratch1, scratch1, 1);
		for (int i = 0; i < Lp*Lp; i++) {
			scratch1[i] *= dt;
			scratch1[i] = Math.max(scratch1[i], -n0[i]/2.); // clip to ensure n(t+dt) > n(t)/2
			n0[i] += scratch1[i]; // n0' <- n0 + dt \del^2 (df/n0)
		}

		//
		// Update n1
		//
		
		fft.convolve(n1, scratch1, potential11);
		fft.convolve(n0, scratch2, potential01);
		for (int i = 0; i < Lp*Lp; i++) {
			scratch1[i] = r*n1[i] + scratch1[i] + cpl*scratch2[i] + n1[i]*n1[i]*n1[i] - trim*1/n1[i]; 
		}
		laplacian(scratch1, scratch1, 1);
		for (int i = 0; i < Lp*Lp; i++) {
			scratch1[i] *= dt;
			scratch1[i] = Math.max(scratch1[i], -n1[i]/2.); // clip to ensure n(t+dt) > n(t)/2
			n1[i] += scratch1[i]; // n1' <- n1 + dt \del^2 (df/n1)
		}

		
		t += dt;
	}
	
	public double energy() {
		return 0;
	}

	private double noise() {
		return noise*random.nextGaussian();
	}
}



public class BinaryPFC2dApp extends Simulation {
	Grid grid1 = new Grid("Grid");
	Grid grid2 = new Grid("Grid");


	BinaryPFC2d sim;
	
	public static void main(String[] args) {
		new Control(new BinaryPFC2dApp(), "Binary Phase field crystal model");
	}

	public void load(Control c) {
		c.frame(grid1);
		c.frame(grid2);
		params.add("L", 20);
		params.add("dx", 0.2);
		params.add("rho0", new DoubleValue(0.5, -1, 1));
		params.add("rho1", new DoubleValue(0.5, -1, 1));
		params.addm("r", new DoubleValue(0, -2, 0.5).withSlider());
		params.addm("a0", new DoubleValue(1, 0.5, 2).withSlider());
		params.addm("a1", new DoubleValue(1, 0.5, 2).withSlider());
		params.addm("dt", 0.5);
		params.addm("Noise", new DoubleValue(0.01, 0, 1).withSlider());
		params.add("Time");
		params.add("Energy");
	}
	
	public void animate() {
		sim.readParams(params);
		
		grid1.setAutoScale();
		grid1.setDrawRange(true);
		grid2.setAutoScale();
		grid2.setDrawRange(true);
		
		int Lp = sim.Lp;
		grid1.registerData(Lp, Lp, sim.n0);
		grid2.registerData(Lp, Lp, sim.n1);

		params.set("Time", format(sim.t));
		params.set("Energy", format(sim.energy()));
	}
	
	public void clear() {
		grid1.clear();
		grid2.clear();
	}
	
	public void run() {
		sim = new BinaryPFC2d(params);
		// sim.randomize();
		Job.animate();
		
		while (true) {
			for (int i = 0; i < 5; i++)
				sim.simulate();
			Job.animate();
		}
	}

}
