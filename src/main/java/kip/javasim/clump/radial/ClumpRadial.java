package kip.javasim.clump.radial;

import static java.lang.Math.PI;
import static java.lang.Math.acos;
import static java.lang.Math.log;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.sqrt;
import static scikit.numerics.Math2.sqr;
import kip.javasim.clump.dim2.AbstractClump2D;
import kip.javasim.clump.dim3.AbstractClump3D;
import scikit.jobs.params.Parameters;
import scikit.util.DoubleArray;

//import static kip.util.MathPlus.hypot;
//import kip.clump.dim3.AbstractClump3D;
//import kip.util.MathPlus;
//import scikit.dataset.Accumulator;
//import scikit.numerics.fft.util.FFT3D;
//import scikit.numerics.fn.Function3D;
//import scikit.util.Commands;

public class ClumpRadial {
	double DENSITY = 1;
	int dim;
	double dx, L, R, T, t, dt;
	int Lp;
	double[] phi, phibar, del_phi;
	double[] freeEnergy, sqr_dF_dphi;
	double rms_dF_dphi, freeEnergyDensity;
	double[][] shell;
	
	
	public ClumpRadial(Parameters params) {
		dim = params.iget("Dimension");
		dx = params.fget("dx");
		L = params.fget("L");
		R = params.fget("R");
		T = params.fget("T");
		dt = params.fget("dt");
		
		Lp = (int)(L/dx);
		dx = L/Lp;
		params.set("dx", dx);
		
		phi = new double[Lp];
		phibar = new double[Lp];
		del_phi = new double[Lp];
		freeEnergy = new double[Lp];
		sqr_dF_dphi = new double[Lp];
		
		for (int i = 0; i < Lp; i++) {
			// phi[i] = (i*dx < R/2.) ? 1 : 0;
			phi[i] = 0.1/(1+sqr(i*dx/R)) + DENSITY;
		}
		t = 0;
		
		buildConvolutionMatrix(R);
	}
	
	public void readParams(Parameters params) {
		T = params.fget("T");
		dt = params.fget("dt");
	}
	
	public void simulate() {
		convolveWithRange(phi, phibar, R);

		for (int ap = 1; ap < Lp; ap++) {
			del_phi[ap] = - dt*(phibar[ap]+T*log(phi[ap]));
		}
		double mu = rmean(del_phi)-(DENSITY-rmean(phi));
		for (int ap = 1; ap < Lp; ap++) {
			// clip del_phi to ensure phi(t+dt) > phi(t)/2
			del_phi[ap] = max(del_phi[ap]-mu, -phi[ap]/2.);
		}
		
		for (int ap = 1; ap < Lp; ap++) {
			sqr_dF_dphi[ap] = sqr(del_phi[ap] / dt);
			freeEnergy[ap] = 0.5*phi[ap]*phibar[ap]+T*phi[ap]*log(phi[ap]);
			phi[ap] += del_phi[ap];
		}
		phi[0] = phi[1]; // for aesthetics
		
		rms_dF_dphi = sqrt(rmean(sqr_dF_dphi));
		freeEnergyDensity = rmean(freeEnergy) - 0.5;
		t += dt;
	}
	
	public void scaleField(double scale) {
		// phi will not be scaled above PHI_UB or below PHI_LB
		double PHI_UB = 5;
		double PHI_LB = 0.01;
		double s1 = (PHI_UB-DENSITY)/(DoubleArray.max(phi)-DENSITY+1e-10);
		double s2 = (PHI_LB-DENSITY)/(DoubleArray.min(phi)-DENSITY-1e-10);
		boolean rescaleClipped = scale > min(s1,s2);
		if (rescaleClipped)
			scale = min(s1,s2);
		for (int ap = 0; ap < Lp; ap++) {
			phi[ap] = scale*(phi[ap] - DENSITY) + DENSITY;
		}
	}
	
	public double rmean(double[] src) {
		double acc = 0;
		double vol = 0;
		for (int ap = 1; ap < Lp; ap++) {
			double shell = (dim == 2) ? 2*PI*ap : 4*PI*ap*ap;
			acc += shell*src[ap];
			vol += shell;
		}
		// 2D: vol -> Pi Lp^2 as dx -> 0
		// 3D: vol -> 4/3 Pi Lp^3
		return acc / vol;
	}
	
	public double rvariance(double[] src) {
		double acc = 0;
		double vol = 0;
		for (int ap = 1; ap < Lp; ap++) {
			double shell = (dim == 2) ? 2*PI*ap : 4*PI*ap*ap;
			acc += shell*sqr(src[ap]);
			vol += shell;
		}
		// 2D: vol -> Pi Lp^2 as dx -> 0
		// 3D: vol -> 4/3 Pi Lp^3
		return acc/vol - sqr(rmean(src));
	}
	
	public double reducedTemperature() {
		double T_sp = (dim == 2) ? AbstractClump2D.T_SP : AbstractClump3D.T_SP;
		return (T-T_sp)/T_sp;
	}
	
	public double KR_SP() {
		return (dim == 2) ? AbstractClump2D.KR_SP : AbstractClump3D.KR_SP;
	}
	
	
	void buildConvolutionMatrix(double R) {
		shell = new double[Lp][];
		int Rp = (int)(R/dx);
		for (int ap = 1; ap < Lp; ap++) {
			int b0 = max(1, ap-Rp);
			int bf = ap+Rp;
			shell[ap] = new double[bf-b0];
			for (int bp = b0; bp < bf; bp++) {
				double a = ap*dx;
				double b = bp*dx;
				double cos_theta = min(max((a*a+b*b-R*R)/(2*a*b), -1), 1);
				if (dim == 2) {
					double theta = acos(cos_theta);
					shell[ap][bp-b0] = dx * 2*theta*b;
				}
				else {
					shell[ap][bp-b0] = dx * (1-cos_theta)*2*PI*b*b;
				}
			}
		}
	}
	
	public void convolveWithRange(double[] src, double[] dst, double R) {
		int Rp = (int)(R/dx);
		for (int ap = 1; ap < Lp; ap++) {
			dst[ap] = 0;
			double vol = 0;
			int b0 = max(1, ap-Rp);
			int bf = ap+Rp;
			for (int bp = b0; bp < bf; bp++) {
				dst[ap] += shell[ap][bp-b0] * (bp >= Lp ? DENSITY : src[bp]);
				vol += shell[ap][bp-b0];
			}
			// 2D: vol -> Pi R^2 as dx -> 0
			// 3D: vol -> 4/3 Pi R^3
			dst[ap] /= vol;
		}
		dst[0] = dst[1]; // for aesthetics
	}

//	public void convolveWithRange2(double[] src, double[] dst, double R) {
//		int Rp = (int)(R/dx);
//		for (int ap = 1; ap < Lp; ap++) {
//			dst[ap] = 0;
//			double vol = 0;
//			for (int bp = max(1, ap-Rp); bp < ap+Rp; bp++) {
//				double a = ap*dx;
//				double b = bp*dx;
//				double cos_theta = min(max((a*a+b*b-R*R)/(2*a*b), -1), 1);
//				double shell;
//				if (dim == 2) {
//					double theta = acos(cos_theta);
//					shell = dx * 2*theta*b;
//				}
//				else {
//					shell = dx * (1-cos_theta)*2*PI*b*b;
//				}
//				dst[ap] += shell * (bp >= Lp ? DENSITY : src[bp]);
//				vol += shell;
//			}
//			// 2D: vol -> Pi R^2 as dx -> 0
//			// 3D: vol -> 4/3 Pi R^3
//			dst[ap] /= vol;
//		}
//		dst[0] = dst[1]; // for aesthetics
//	}
	
	
//	int fLp = 128;
//	double[] scratch = new double[fLp*fLp*fLp];
//	FFT3D fft = new FFT3D(fLp, fLp, fLp);
//	
//	public void convolveWithRangeSlow(double[] src, double[] dst, final double R) {
//		double fL = L*2.5;
//		double fdx = fL/fLp;
//		
//		for (int zp = 0; zp < fLp; zp++) {
//			for (int yp = 0; yp < fLp; yp++) {
//				for (int xp = 0; xp < fLp; xp++) {
//					double x = fdx*(xp-fLp/2);
//					double y = fdx*(yp-fLp/2);
//					double z = fdx*(zp-fLp/2);
//					double r = MathPlus.hypot(x, y, z);
//					int ip = (int)(r/dx);
//					double phi = (ip >= Lp) ? DENSITY : src[ip];
//					scratch[zp*fLp*fLp + yp*fLp + xp] = phi;
//				}
//			}
//		}
//		
//		fft.setLengths(fL,fL,fL);
//		fft.convolve(scratch, scratch, new Function3D() {
//			public double eval(double k1, double k2, double k3) {
//				return AbstractClump3D.potential(hypot(k1*R,k2*R,k3*R));
//			}
//		});
//		
//		// Commands.grid3d(fLp,fLp,fLp,scratch);
//
//		Accumulator acc = new Accumulator(1*dx);
//		acc.setAveraging(true);
//		for (int zp = 0; zp < fLp; zp++) {
//			for (int yp = 0; yp < fLp; yp++) {
//				for (int xp = 0; xp < fLp; xp++) {
//					double x = fdx*(xp-fLp/2);
//					double y = fdx*(yp-fLp/2);
//					double z = fdx*(zp-fLp/2);
//					double r = MathPlus.hypot(x, y, z);
//					double phi = scratch[zp*fLp*fLp + yp*fLp + xp];
//					acc.accum(r, phi);
//				}
//			}
//		}
//		
//		for (int ap = 0; ap < Lp; ap++) {
//			double v = acc.eval(ap*dx);
//			dst[ap] = Double.isNaN(v) ? dst[ap-1] : v;
//		}
//		
//		// Commands.plot(dst);		
//	}
	
}
