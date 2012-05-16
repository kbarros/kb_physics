package kip.javasim.ising.dim1;

import scikit.jobs.params.Parameters;
import scikit.numerics.fft.managed.ComplexDoubleFFT;
import scikit.numerics.fft.managed.ComplexDoubleFFT_Mixed;
import static java.lang.Math.*;
import static scikit.numerics.Math2.*;


public class FieldIsing extends AbstractIsing {
	public double[] phi, phibar, phi2, phibar2;
	ComplexDoubleFFT fft;	// Object to perform transforms
	public double[] fftScratch;
	int Lp;
	
	public FieldIsing(Parameters params) {
		initialize(params);
	}
	
	
    public FieldIsing clone() {
		FieldIsing c = (FieldIsing)super.clone();
		c.phi = (double[])phi.clone();
		return c;
    }
	
	
	// reset time, set random number seed, initialize fields to down
	public void initialize(Parameters params) {
		super.initialize(params);
		Lp = N/dx;
		phi = new double[Lp];
		phibar = new double[Lp];
		phi2 = new double[Lp];
		phibar2 = new double[Lp];
		
		fftScratch = new double[2*Lp];
		fft = new ComplexDoubleFFT_Mixed(Lp);		
	}
	
	
	public void randomizeField(double m) {
		for (int i = 0; i < Lp; i++)
			phi[i] = m + random.nextGaussian()*sqrt((1-m*m)/dx);
	}
	
	public void setField(double m) {
		for (int i = 0; i < Lp; i++)
			phi[i] = m;
	}
	
	public double magnetization() {
		double sum = 0;
		for (int i = 0; i < Lp; i++)
			sum += phi[i];
		return sum / Lp;
	}
	
	
	public double fieldElement(int i) {
		return phi[i];
	}
	
	
	private void convolveWithRange(double[] src, double[] dest) {
		for (int i = 0; i < Lp; i++) {
			fftScratch[2*i] = src[i];
			fftScratch[2*i+1] = 0;
		}
		
		fft.transform(fftScratch);
		for (int x = -Lp/2; x < Lp/2; x++) {
				double kR = (2*PI*x/N) * R;
				int i = (x+Lp)%Lp;
				double J = (kR == 0 ? 1 : sin(kR)/kR);
				fftScratch[2*i] *= J;
				fftScratch[2*i+1] *= J;
		}
		fft.backtransform(fftScratch);
		
		for (int i = 0; i < Lp; i++) {
			dest[i] = fftScratch[2*i] / Lp;
		}		
	}
	
	
//	private double bar(double[] field, int i) {
//		int nblocks = (int)round(2.0*R/dx);
//		
//		double acc = field[i];
//		for (int j = 1; j <= nblocks/2; j++) {
//			acc += field[(i+j)%Lp];
//			acc += field[(i-j+Lp)%Lp];
//		}
//		return acc / (1 + 2*(nblocks/2));
//	}
	
	private double drift(double phi, double phibar) {
		phibar = (J*phibar+h)/T;
		double tphi = tanh(phibar);
		switch (dynamics) {
		case GLAUBER:
			return tphi - phi;
		case METROPOLIS:
			return 2*exp(-abs(phibar)) * (sinh(phibar) - phi*cosh(phibar));
		default:
			return Double.NaN;
		}
	}
	
	public double noise(double phi, double phibar) {
		double tphi = tanh((J*phibar+h)/T);
		switch (dynamics) {
		case GLAUBER:
			return sqrt(2 - tphi*tphi - phi*phi);
		case METROPOLIS:
			return 2; // ???
		default:
			return Double.NaN;
		}
	}
	
	protected void _step() {
		// break update into "steps" number of substeps. each one with time interval dt_
		double dx_ = dx;
		double CONST = 100000;
		int steps = (int) max(dt*CONST/(dx_*dx_), 1);
		double dt_ = dt / steps; // = dx_/k
		
		// full heun scheme
		for (int cnt = 0; cnt < steps; cnt++) {
			convolveWithRange(phi, phibar);
			
			// get euler predictor
			for (int i = 0; i < N/dx; i++) {
				if (sqr(phi[i]) > 1)
					throw new ArithmeticException();
				double U = drift(phi[i], phibar[i]);
				double V = noise(phi[i], phibar[i]);
				double eta = random.nextGaussian();
				phi2[i] = phi[i] + dt_*U + sqrt(dt_/dx_)*eta*V;
			}
			
			convolveWithRange(phi2, phibar2);
			
			// take step based on predictor
			for (int i = 0; i < N/dx; i++) {
				double U1 = drift(phi[i],  phibar[i]);
				double U2 = drift(phi2[i], phibar2[i]);
				double V1 = noise(phi[i],  phibar[i]);
				double V2 = noise(phi2[i], phibar2[i]);
				double eta = random.nextGaussian();
				phi[i] += dt_*(U1+U2)/2 + sqrt(dt_/dx_)*eta*(V1+V2)/2;
			}
		}
	}
}