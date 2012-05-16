package kip.javasim.langevin.apps;

import kip.javasim.Random;
import scikit.jobs.params.Parameters;
import scikit.numerics.fft.FFT3D;
import scikit.numerics.fft.FFT3DManaged;
import scikit.numerics.fn.Function3D;


class TDGL3d {
	Random random;
	double L, dx, dt, t;
	int Lp;
	
	// constants from Cheng & Rutenberg
	double a1 = 3;
	double a2 = 0;
	
	// approximate energy at a single interface per unit area
	double r = 1;
	double surfaceEnergyDensity = 1.0; // 0.94;
	
	double[] phi, scratch1, scratch2;
	FFT3D fft;
	String kernel;
	
	public TDGL3d(Parameters params) {
		random = new Random(params.iget("Random seed", 0));
		
		L = params.fget("L");
		dx = params.fget("dx");
		dt = params.fget("dt");

		Lp = (int) (L / dx);
		dx = L / Lp;
		params.set("dx", dx);
		
		phi = new double[Lp*Lp*Lp];
		scratch1 = new double[Lp*Lp*Lp];
		scratch2 = new double[Lp*Lp*Lp];
		
		fft = new FFT3DManaged(Lp, Lp, Lp);
		fft.setLengths(L, L, L);
		
		randomize();
	}
	
	public void randomize() {
		for (int i = 0; i < Lp*Lp*Lp; i++)
			phi[i] = 0.1*random.nextGaussian();
		t = 0;
	}
	
	public void randomizeAndShift() {
		randomize();
		double sum = 0;
		for (int i = 0; i < Lp*Lp*Lp; i++)
			sum += phi[i];
		for (int i = 0; i < Lp*Lp*Lp; i++)
			phi[i] -= sum / (Lp*Lp*Lp);
		
		sum = 0;
		for (int i = 0; i < Lp*Lp*Lp; i++)
			sum += phi[i];
		if (sum > 1e-8 || sum < -1e-8) {
			System.out.println("Error in shifting!");
		}

		t = 0;
	}

	public void randomizeIsing() {
		int Lp3 = Lp*Lp*Lp; 
		int[] indices = new int[Lp3];
		for (int i = 0; i < Lp3; i++)
			indices[i] = i;
		// randomize indices
		for (int i = 0; i < Lp3; i++) {
			int r = i+random.nextInt(Lp3-i);
			// swap indices[i] and indices[r]
			int that = indices[r];
			indices[r] = indices[i];
			indices[i] = that;
		}
		for (int i = 0; i < Lp3; i++) {
			phi[indices[i]] = (i < Lp3/2) ? -1 : +1;
		}

		double sum = 0;
		for (int i = 0; i < Lp3; i++)
			sum += phi[i];
		if (sum != 0.0) {
			System.out.println("Error in randomizing!");
		}
		
		t = 0;
	}
	
	public void schwartzPSurface() {
		// initial condition designed to evolve into schwarz p surface
		for (int x = 0; x < Lp; x++) {
			for (int y = 0; y < Lp; y++) {
				for (int z = 0; z < Lp; z++) {
					int ix = (x < Lp/2) ? 0 : 1;
					int iy = (y < Lp/2) ? 0 : 1;
					int iz = (z < Lp/2) ? 0 : 1;
					phi[z*Lp*Lp + y*Lp + x] += (ix + iy + iz <= 1) ? 1 : -1;
				}
			}
		}
	}
	
	public void readParams(Parameters params) {
		dt = params.fget("dt");
	}
		
	public double r2k2(double kx, double ky, double kz) {
		double k2 = kx*kx + ky*ky + kz*kz;
		double c = 27*(kx*kx*ky*ky*kz*kz)/(k2*k2*k2);
		double a = 0.90;
		return k2==0 ? 0 : r*r*((1-a)+a*c)*k2;
	}

	public double netMagnetization() {
		double sum = 0;
		for (int i = 0; i < phi.length; i++)
			sum += phi[i];
		return sum;
	}
	
	public void simulate() {
		System.out.println(netMagnetization());
		
		// scratch1 stores term proportional to phi
		fft.convolve(phi, scratch1, new Function3D() {
			public double eval(double kx, double ky, double kz) {
				//double k2 = kx*kx + ky*ky + kz*kz;
				double k2 = r2k2(kx, ky, kz);
				return (1 + a1*dt + a2*dt*(-k2)) / (1 + (a1-1)*dt + (a2-1)*dt*(-k2));
			}
		});
		
		// scratch2 stores term proportional to phi^3
		for (int i = 0; i < Lp*Lp*Lp; i++) {
			scratch2[i] = phi[i]*phi[i]*phi[i];
		}
		fft.convolve(scratch2, scratch2, new Function3D() {
			public double eval(double kx, double ky, double kz) {
				//double k2 = kx*kx + ky*ky + kz*kz;
				double k2 = r2k2(kx, ky, kz);
				return dt / (1 + (a1-1)*dt + (a2-1)*dt*(-k2));
			}
		});

		for (int i = 0; i < Lp*Lp*Lp; i++) {
			phi[i] = scratch1[i] - scratch2[i];
		}
		
		t += dt;
	}
	
	// F = \int dx [ - phi (1 + \del^2) phi / 2 + phi^4 / 4 ]  
	public double freeEnergy() {
		// scratch1 stores (1 + \del^2) phi
		fft.convolve(phi, scratch1, new Function3D() {
			public double eval(double kx, double ky, double kz) {
				//double k2 = kx*kx + ky*ky + kz*kz;
				double k2 = r2k2(kx, ky, kz);
				return 1 - k2;
			}
		});
		double ret = 0;
		for (int i = 0; i < Lp*Lp*Lp; i++) {			
			ret += - phi[i]*scratch1[i]/2;      // - phi (1 + \del^2) phi / 2 
			ret += phi[i]*phi[i]*phi[i]*phi[i] / 4;	// phi^4 / 4
		}
		
		// shift energy so that it is zero in the ground state
		ret += Lp*Lp*Lp/4; 
		
		return ret*dx*dx*dx;
	}

}