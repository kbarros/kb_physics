package kip.javasim.ising.dim2;

import static java.lang.Math.rint;
import static java.lang.Math.sqrt;
import static scikit.numerics.Math2.hypot;
import static scikit.numerics.Math2.sqr;
import kip.javasim.Random;
import scikit.dataset.Accumulator;
import scikit.jobs.params.Parameters;
import scikit.numerics.fft.FFT2D;
import scikit.numerics.fn.Function2D;


class PhiFourth2DConsts {
	static public Function2D potential1 = new Function2D() {
		public double eval(double kx, double ky) {
			return - (kx*kx + ky*ky);
		}
	};
	
	static public Function2D potential2 = new Function2D() {
		public double eval(double kx, double ky) {
			double theta = Math.atan2(ky, kx);
			return - (kx*kx + ky*ky)*(1+Math.cos(4*theta)/2.);
		}
	};
}

public class PhiFourth2D {

	public double L, R, T, h, dx;
	Random random = new Random();
	
	int Lp;
	double t;
	double[] phi, laplace_phi, del_phi, potential;
	FFT2D fft;
	boolean noiselessDynamics = false;
	
	public double dt;
	public double rms_dF_dphi;
	public double freeEnergyDensity;
	public double anisotropy;
	
	
	public PhiFourth2D(Parameters params) {
		random.setSeed(params.iget("Random seed", 0));
		
		R = params.fget("R");
		L = R*params.fget("L/R");
		T = params.fget("T");
		h = params.fget("h");
		dx = R*params.fget("dx/R");
		dt = params.fget("dt");
		anisotropy = params.fget("anisotropy");
		
		noiselessDynamics = params.sget("Noise").equals("No");
		Lp = Integer.highestOneBit((int)rint(L/dx));
		dx = L / Lp;
		params.set("dx/R", dx/R);

		phi = new double[Lp*Lp];
		laplace_phi = new double[Lp*Lp];
		del_phi = new double[Lp*Lp];
		
		fft = new FFT2D(Lp, Lp);
		fft.setLengths(L, L);
		
		potential = fft.buildFourierArray(PhiFourth2DConsts.potential2);
		
		t = 0;
		for (int i = 0; i < Lp*Lp; i++)
			phi[i] = 0;
	}
	
	public double[] phi() {
		return phi;
	}
	
	public void readParams(Parameters params) {
		noiselessDynamics = params.sget("Noise").equals("No");
		T = params.fget("T");
		h = params.fget("h");
		dt = params.fget("dt");
		anisotropy = params.fget("anisotropy");
	}
	
	
	public void stripes() {
		for (int x = 0; x < Lp; x++) {
			for (int y = 0; y < Lp; y++) {
				int i = y*Lp + x;
				phi[i] = Math.sin(0.5*(x+32)); // hypot(x*dx, y*dx); 
			}
		}		
	}
	
	public void randomize() {
		for (int i = 0; i < Lp*Lp; i++) {
			phi[i] = sqrt(1/(dx*dx))*noise();
		}
	}
	
	public void useNoiselessDynamics(boolean b) {
		noiselessDynamics = b;
	}	
	
//	private void laplaceOperator(double[] phi, double[] phi_bar) {
//		for (int x = 0; x < Lp; x++) {
//			for (int y = 0; y < Lp; y++) {
//				int xp = (x+1)%Lp;
//				int xm = (x-1+Lp)%Lp;
//				int yp = (y+1)%Lp;
//				int ym = (y-1+Lp)%Lp;
//				phi_bar[y*Lp+x] = (-4*phi[y*Lp+x] + phi[yp*Lp+x] + phi[ym*Lp+x] + 
//									phi[y*Lp+xp] + phi[y*Lp+xm]) / (dx*dx); 
//			}
//		}
//	}
	
	public void simulate() {
//		laplaceOperator(phi, laplace_phi);
//		fft.convolve(phi, laplace_phi, PhiFourth2DConsts.potential2);
		
		Function2D potentialFunc = new Function2D() {
			public double eval(double kx, double ky) {
				double theta = Math.atan2(ky, kx);
				return - (kx*kx + ky*ky)*(1+anisotropy*Math.cos(4*theta)/2.);
			}
		};
		fft.buildFourierArray(potential, potentialFunc);
		fft.convolve(phi, laplace_phi, potential);
		
		for (int i = 0; i < Lp*Lp; i++) {
			del_phi[i] = - dt*(-R*R*laplace_phi[i] + (T+sqr(phi[i]))*phi[i] - h);
			if (!noiselessDynamics)
				del_phi[i] += sqrt(dt/(dx*dx))*noise();
		}
		
		freeEnergyDensity = 0;
		for (int i = 0; i < Lp*Lp; i++) {
			freeEnergyDensity += -phi[i]*(R*R*laplace_phi[i])/2 + sqr(T+sqr(phi[i]))/4 - h*phi[i];
			phi[i] += del_phi[i];
		}
		freeEnergyDensity /= (Lp*Lp);
		t += dt;
	}

	public Accumulator newStructureAccumulator(double binWidth) {
		return new Accumulator(binWidth);
	}
	
	public void accumulateStructure(final Accumulator sf) {
		fft.transform(phi, new FFT2D.MapFn() {
			public void apply(double k1, double k2, double re, double im) {
				double kmag = hypot(k1, k2);
				if (kmag > 0 && kmag <= 4)
					sf.accum(kmag, (re*re+im*im)/(L*L));
			}
		});
	}	
	
	public int numColumns() {
		return Lp;
	}
	
	public double time() {
		return t;
	}
	
	private double noise() {
		return random.nextGaussian();
	}

}
