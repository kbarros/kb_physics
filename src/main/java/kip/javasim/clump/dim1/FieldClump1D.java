package kip.javasim.clump.dim1;

import static java.lang.Math.floor;
import static java.lang.Math.log;
import static java.lang.Math.max;
import static java.lang.Math.rint;
import static java.lang.Math.sqrt;
import static scikit.numerics.Math2.sqr;
import scikit.dataset.Accumulator;
import scikit.jobs.params.Parameters;
import scikit.numerics.fft.FFT1D;
import scikit.numerics.fn.Function1D;
import scikit.util.DoubleArray;

public class FieldClump1D extends AbstractClump1D {
	int Lp;
	double dt, t;
	double[] phi, phi_bar, del_phi;
	FFT1D fft;
	
	boolean noiselessDynamics = false;
	public double rms_dF_dphi;
	public double freeEnergyDensity;
	
	
	public FieldClump1D(Parameters params) {
		random.setSeed(params.iget("Random seed", 0));
		
		R = params.fget("R");
		L = params.fget("L");
		T = params.fget("T");
		dx = params.fget("dx");
		dt = params.fget("dt");
		
		Lp = Integer.highestOneBit((int)rint(L/dx));
		dx = L / Lp;
		params.set("dx", dx);
		allocate();
		
		t = 0;
		for (int i = 0; i < Lp; i++)
			phi[i] = DENSITY;
	}
	
	public void readParams(Parameters params) {
		T = params.fget("T");
		dt = params.fget("dt");
	}
	
	
	public void initializeFieldWithSeed(String type) {
  		for (int i = 0; i < Lp; i++) {
			double x = dx*(i - Lp/2);
			double field = random.nextGaussian();
			double mag = 0.1 / (1+sqr(x/R));
			phi[i] = DENSITY*(1+mag*field);
		}
	}
	
	public void useNoiselessDynamics(boolean b) {
		noiselessDynamics = b;
	}
	
	public void simulate() {
		fft.convolve(phi, phi_bar, new Function1D() {
			public double eval(double k) {
				return potential(k*R);
			}
		});
		
		for (int i = 0; i < Lp; i++) {
			del_phi[i] = - dt*(phi_bar[i]+T*log(phi[i])) + sqrt(dt*2*T/dx)*noise();
		}
		double mu = DoubleArray.mean(del_phi)-(DENSITY-DoubleArray.mean(phi));
		for (int i = 0; i < Lp; i++) {
			// clip del_phi to ensure phi(t+dt) > phi(t)/2
			del_phi[i] = max(del_phi[i]-mu, -phi[i]/2.);
		}
		
		rms_dF_dphi = 0;
		freeEnergyDensity = 0;
		for (int i = 0; i < Lp; i++) {
			rms_dF_dphi += sqr(del_phi[i] / dt);
			freeEnergyDensity += 0.5*phi[i]*phi_bar[i]+T*phi[i]*log(phi[i]);
			phi[i] += del_phi[i];
		}
		rms_dF_dphi = sqrt(rms_dF_dphi/Lp);
		freeEnergyDensity /= Lp;
		freeEnergyDensity -= 0.5;
		t += dt;
	}
	
	public double dFdensity_dR() {
		double[] dphibar_dR = phi_bar;
		fft.convolve(phi, dphibar_dR, new Function1D() {
			public double eval(double k) {
				return dpotential_dkR(k*R)*k;
			}
		});
		double ret = 0;
		for (int i = 0; i < Lp; i++) {
			ret += 0.5*phi[i]*dphibar_dR[i];
		}
		return ret / Lp;
	}
	
	public Accumulator newStructureAccumulator(double binWidth) {
		// round binwidth down so that it divides KR_SP without remainder.
		binWidth = KR_SP / floor(KR_SP/binWidth);
		return new Accumulator(binWidth);
	}
	
	
	public void accumulateStructure(final Accumulator sf) {
		fft.transform(phi, new FFT1D.MapFn() {
			public void apply(double k, double re, double im) {
				double kR = k*R;
				if (kR > 0 && kR <= 4*KR_SP)
					sf.accum(kR, (re*re+im*im)/L);
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
		phi = new double[Lp];
		phi_bar = new double[Lp];
		del_phi = new double[Lp];
		fft = new FFT1D(Lp);
		fft.setLength(L);
	}
	
	private double noise() {
		return noiselessDynamics ? 0 : random.nextGaussian();
	}
}