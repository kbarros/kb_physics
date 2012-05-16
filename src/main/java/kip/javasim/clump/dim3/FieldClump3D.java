package kip.javasim.clump.dim3;

import static java.lang.Math.cos;
import static java.lang.Math.floor;
import static java.lang.Math.log;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.rint;
import static java.lang.Math.sqrt;
import static scikit.numerics.Math2.hypot;
import static scikit.numerics.Math2.sqr;
import static scikit.util.DoubleArray.mean;
import static scikit.util.DoubleArray.variance;
import scikit.dataset.Accumulator;
import scikit.jobs.params.Parameters;
import scikit.numerics.fft.FFT3D;
import scikit.numerics.fn.Function3D;
import scikit.util.DoubleArray;

public class FieldClump3D extends AbstractClump3D {
	// phi will not be scaled above PHI_UB or below PHI_LB
	double PHI_UB = 19*DENSITY;
	double PHI_LB = 0.001*DENSITY;
	
	int Lp;
	double t;
	double[] phi, phi_bar, del_phi;
	FFT3D fft;
	boolean noiselessDynamics = false;
	
	public double packingFraction = 0;
	public double dt;
	public double Rx, Ry, Rz;
	public boolean rescaleClipped = false; // indicates saddle point invalid
	public double rms_dF_dphi;
	public double freeEnergyDensity;
	
	
	public FieldClump3D(Parameters params) {
		random.setSeed(params.iget("Random seed", 0));
		
		Rx = Ry = Rz = params.fget("R");
		L = params.fget("L");
		T = params.fget("T");
		dx = params.fget("dx");
		dt = params.fget("dt");
		
		Lp = Integer.highestOneBit((int)rint(L/dx));
		dx = L / Lp;
		params.set("dx", dx);
		allocate();
		
		t = 0;
		for (int i = 0; i < Lp*Lp*Lp; i++)
			phi[i] = DENSITY;
	}
	
	public void halveResolution() {
		int old_Lp = Lp;
		double[] old_phi = phi; 
		Lp /= 2;
		dx *= 2.0;
		allocate();
		for (int z = 0; z < Lp; z++) {
			for (int y = 0; y < Lp; y++) {
				for (int x = 0; x < Lp; x++) {
					phi[z*Lp*Lp + y*Lp + x] = old_phi[2*z*old_Lp*old_Lp + 2*y*old_Lp + 2*x];
				}
			}
		}
	}
	
	public void doubleResolution() {
		int old_Lp = Lp;
		double[] old_phi = phi; 
		Lp *= 2;
		dx /= 2.0;
		allocate();
		for (int z = 0; z < Lp; z++) {
			for (int y = 0; y < Lp; y++) {
				for (int x = 0; x < Lp; x++) {
					phi[z*Lp*Lp + y*Lp + x] = old_phi[(z/2)*old_Lp*old_Lp + (y/2)*old_Lp + (x/2)];
				}
			}
		}
	}
	
	public void duplicateAndTile() {
		int Lp2 = Lp;
		double[] old_phi = phi;
		L *= 2;
		Lp *= 2;
		allocate();
		for (int z = 0; z < Lp; z++) {
			for (int y = 0; y < Lp; y++) {
				for (int x = 0; x < Lp; x++) {
					double p = old_phi[(z%Lp2)*Lp2*Lp2 + (y%Lp2)*Lp2 + (x%Lp2)];
					if (true) {
						double r = dx * hypot((z-Lp/2), (y-Lp/2), (x-Lp/2)) / Rx;
						p = (p - DENSITY) / (1 + 0.1*r*r) + DENSITY;
					}
					phi[z*Lp*Lp + y*Lp + x] = p;
				}
			}
		}
	}
	
	public void duplicateAndEmbed() {
		int Lp2 = Lp;
		double[] old_phi = phi;
		L *= 2;
		Lp *= 2;
		allocate();
		DoubleArray.set(phi, DENSITY);
		for (int z = 0; z < Lp2; z++) {
			for (int y = 0; y < Lp2; y++) {
				for (int x = 0; x < Lp2; x++) {
					double p = old_phi[z*Lp2*Lp2 + y*Lp2 + x];
					phi[(z+Lp2/2)*Lp*Lp + (y+Lp2/2)*Lp + (x+Lp2/2)] = p;
				}
			}
		}
	}
	
	public void readParams(Parameters params) {
		T = params.fget("T");
		dt = params.fget("dt");
	}
	
	
	public void initializeFieldWithSeed(String type) {
  		for (int i = 0; i < Lp*Lp*Lp; i++) {
			double R = Rx;
			double x = dx*(i%Lp - Lp/2);
			double y = dx*((i%(Lp*Lp))/Lp - Lp/2);
			double z = dx*((i/(Lp*Lp)) - Lp/2);
			double field = 0;
			double k = KR_SP/R;
			if (type.equals("BCC")) {
				// BCC (reciprocal lattice is FCC)
				field += cos(k * ( x + z) / sqrt(2));
				field += cos(k * (-x + z) / sqrt(2));
				field += cos(k * ( y + z) / sqrt(2));
				field += cos(k * (-y + z) / sqrt(2));
			}
			else if (type.equals("Triangle")) {
				double rad = 0.5*R;
				double sigma = 0.1*R;
				double x0 = rad, y0 = 0, z0 = 0;
				field += Math.exp((-sqr(x-x0)-sqr(y-y0)-sqr(z-z0))/(2*sqr(sigma)));
				x0 = rad*cos(2*Math.PI/3);
				y0 = rad*Math.sin(2*Math.PI/3);
				field += Math.exp((-sqr(x-x0)-sqr(y-y0)-sqr(z-z0))/(2*sqr(sigma)));
				x0 = rad*cos(4*Math.PI/3);
				y0 = rad*Math.sin(4*Math.PI/3);
				field += Math.exp((-sqr(x-x0)-sqr(y-y0)-sqr(z-z0))/(2*sqr(sigma)));
			}
			else if (type.equals("Noise")) {
				// random
				field = random.nextGaussian();
			}
			double r = sqrt(x*x+y*y+z*z);
			double mag = 0.2 / (1+sqr(r/R));
			phi[i] = DENSITY*(1+mag*field);
		}
	}
	
	public void useNoiselessDynamics(boolean b) {
		noiselessDynamics = b;
	}
	
	public double phiVariance() {
		return variance(phi);
	}
	
	public void scaleField(double scale) {
		double s1 = (PHI_UB-DENSITY)/(DoubleArray.max(phi)-DENSITY+1e-10);
		double s2 = (PHI_LB-DENSITY)/(DoubleArray.min(phi)-DENSITY-1e-10);
		rescaleClipped = scale > min(s1,s2);
		if (rescaleClipped)
			scale = min(s1,s2);
		for (int i = 0; i < Lp*Lp*Lp; i++) {
			phi[i] = (phi[i]-DENSITY)*scale + DENSITY;
		}
		
//		rescaleClipped = false;
//		for (int i = 0; i < Lp*Lp*Lp; i++) {
//			phi[i] = (phi[i]-DENSITY)*scale + DENSITY;
//			if (phi[i] < PHI_LB || PHI_UB < phi[i]) {
//				rescaleClipped = true;
//				phi[i] = max(min(phi[i], PHI_UB), PHI_LB);
//			}
//		}
	}
	
	public double entropy(double phi) {
		if (packingFraction > 0) {
			double a = 1/packingFraction;
			return phi*log(phi) + (a-phi)*log(a-phi);
		}
		else {
			return phi*log(phi);
		}
	}
	
	public double dentropy_dphi(double phi) {
		if (packingFraction > 0) {
			double a = 1/packingFraction;
			return log(phi/(a-phi));
		}
		else {
			return log(phi);
		}
	}
	
	public double freeEnergyBackground() {
		return 0.5*sqr(DENSITY) + T*entropy(DENSITY);
	}
	
	public void simulate() {
		fft.convolve(phi, phi_bar, new Function3D() {
			public double eval(double k1, double k2, double k3) {
				return potential(hypot(k1*Rx,k2*Ry,k3*Rz));
			}
		});
		
		for (int i = 0; i < Lp*Lp*Lp; i++) {
			del_phi[i] = - dt*(phi_bar[i]+T*dentropy_dphi(phi[i])) + noise();
		}
		double mu = mean(del_phi)-(DENSITY-mean(phi));
		for (int i = 0; i < Lp*Lp*Lp; i++) {
			// clip del_phi to ensure phi(t+dt) > phi(t)/2
			del_phi[i] = max(del_phi[i]-mu, -phi[i]/2.);
		}
		
		rms_dF_dphi = 0;
		freeEnergyDensity = 0;
		for (int i = 0; i < Lp*Lp*Lp; i++) {
			rms_dF_dphi += sqr(del_phi[i] / dt);
			freeEnergyDensity += 0.5*phi[i]*phi_bar[i]+T*entropy(phi[i]);
			phi[i] += del_phi[i];
		}
		rms_dF_dphi = sqrt(rms_dF_dphi/(Lp*Lp*Lp));
		freeEnergyDensity /= (Lp*Lp*Lp);
		freeEnergyDensity -= freeEnergyBackground();
		t += dt;
	}
	
	public double dFdensity_dRx() {
		double[] dphibar_dR = phi_bar;
		fft.convolve(phi, phi_bar, new Function3D() {
			public double eval(double k1, double k2, double k3) {
				double kR = hypot(k1*Rx, k2*Ry, k3*Rz);
				double dkR_dRx = k1 == 0 ? 0 : (k1*k1*Rx / kR);
				return dpotential_dkR(kR)*dkR_dRx;
			}
		});
		return DoubleArray.dot(phi, dphibar_dR) / (2*Lp*Lp*Lp);
	}
	
	public double dFdensity_dRy() {
		double[] dphibar_dR = phi_bar;
		fft.convolve(phi, phi_bar, new Function3D() {
			public double eval(double k1, double k2, double k3) {
				double kR = hypot(k1*Rx, k2*Ry, k3*Rz);
				double dkR_dRy = k2 == 0 ? 0 : (k2*k2*Ry / kR);
				return dpotential_dkR(kR)*dkR_dRy;
			}
		});
		return DoubleArray.dot(phi, dphibar_dR) / (2*Lp*Lp*Lp);
	}
	
	public double dFdensity_dRz() {
		double[] dphibar_dR = phi_bar;
		fft.convolve(phi, phi_bar, new Function3D() {
			public double eval(double k1, double k2, double k3) {
				double kR = hypot(k1*Rx, k2*Ry, k3*Rz);
				double dkR_dRz = k3 == 0 ? 0 : (k3*k3*Rz / kR);
				return dpotential_dkR(kR)*dkR_dRz;
			}
		});
		return DoubleArray.dot(phi, dphibar_dR) / (2*Lp*Lp*Lp);
	}
	
	public Accumulator newStructureAccumulator(double binWidth) {
		// round binwidth down so that it divides KR_SP without remainder.
		binWidth = KR_SP / floor(KR_SP/binWidth);
		return new Accumulator(binWidth);
	}
	
	public void accumulateStructure(final Accumulator sf) {
		fft.transform(phi, new FFT3D.MapFn() {
			public void apply(double k1, double k2, double k3, double re, double im) {
				double kR = hypot(k1*Rx, k2*Ry, k3*Rz);
				if (kR > 0 && kR <= 4*KR_SP)
					sf.accum(kR, (re*re+im*im)/(L*L*L));
			}
		});
	}
	
	public double[] coarseGrained() {
		return phi;
	}
	
	
	public int numColumns() {
		return Lp;
	}
	
	
	public double time() {
		return t;
	}
	
	private void allocate() {
		phi = new double[Lp*Lp*Lp];
		phi_bar = new double[Lp*Lp*Lp];
		del_phi = new double[Lp*Lp*Lp];
		fft = FFT3D.create(Lp, Lp, Lp);
		fft.setLengths(L, L, L);
	}
	
	private double noise() {
		return noiselessDynamics ? 0 : sqrt(dt*2*T/(dx*dx*dx))*random.nextGaussian();
	}
}