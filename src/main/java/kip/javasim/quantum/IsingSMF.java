package kip.javasim.quantum;

import static java.lang.Math.*;
import static scikit.numerics.Math2.*;
import static scikit.util.DoubleArray.*;

import scikit.dataset.Accumulator;
import scikit.numerics.fft.managed.RealDoubleFFT;
import scikit.numerics.fft.managed.RealDoubleFFT_Radix2;
import scikit.util.Terminal;

class IsingSolution {
	static int N = 1024;
	double J_eff;
	double coef_a;
	double coef_e;
	double J;
	double E0;	// energy per unit volume.
	double[] M_k, K_x;
	
	public IsingSolution(double J_eff) {
		this.J_eff = J_eff;
		coef_a = find_a(J_eff);
		coef_e = coupling_e(J_eff, coef_a);
		J = coef_e*J_eff;
		E0 = E0(J_eff, coef_a, coef_e);
		
		M_k = M_k(J_eff, coef_a);
		K_x = K_x(J_eff, coef_a);
	}
	
	public String toString() {
		return
			"J_eff = " + J_eff +
			"\na = " + coef_a +
			"\ne = " + coef_e +
			"\nJ = " + J +
			"\nE0 = " + E0 +
			"\n";
	}
	
	static double find_a(double J_eff) {
		double a_lo = -0.2;
		double a_hi = 0;
		
		while (a_hi - a_lo > 1e-8) {
			double a = (a_lo + a_hi) / 2.;
			double mm = meanM(J_eff, a);
			if (Double.isNaN(mm) || mm < 0.)
				a_hi = a;
			else
				a_lo = a;
		}
		
		if (Double.isNaN(meanM(J_eff, a_hi))) {
			throw new IllegalArgumentException("Illegal J_eff = " + J_eff);
		}
		return (a_lo + a_hi) / 2.;
	}
	
	static double meanM(double J_eff, double a) {
		return mean(M_k(J_eff, a));
	}
	
	static double[] M_k(double J_eff, double a) {
		double[] m = new double[N];
		for (int k_i = 0; k_i < N; k_i++) {
			double k = 2.*PI*k_i / N;
			m[k_i] = - 1. + sqrt(1 + 4*J_eff*cos(k) - 2.*a);
		}
		return m;
	}
	
	static double[] K_x(double J_eff, double a) {
		double[] M_x = fftReverse(M_k(J_eff, a));
		double[] ret = new double[N];
		for (int i = 0; i < N; i++) {
			ret[i] = atanh(M_x[i]);
		}
		return ret;
	}
	
	static double coupling_e(double J_eff, double a) {
		double[] K = K_x(J_eff, a);
		double ret = 1.;
		for (int i = 0; i < N; i++) {
			ret *= cosh(K[i]);
		}
		return ret;
	}
	
	double E0(double J_eff, double a, double e) {
		double[] M = M_k(J_eff, a);
		
		double tr_M2 = 0;
		for (int i = 0; i < N; i++) {
			tr_M2 += M[i]*M[i] / N;
		}
		return - e * (1. - a - tr_M2 / 2.); 
	}
	
	static double[] fftReverse(double[] a) {
		RealDoubleFFT transf = new RealDoubleFFT_Radix2(N);
		double[] ret = scikit.util.DoubleArray.clone(a);
		transf.transform(ret);
		for (int i = 1; i < N/2; i++) {
			ret[N-i] = ret[i];
		}
		scale(ret, 1./N);
		return ret;
	}
}

public class IsingSMF {
	Terminal term;

	public IsingSolution sol(double J_eff) {
		return new IsingSolution(J_eff);
	}
	
	@SuppressWarnings("unused")
	public Accumulator ground() {
		Accumulator ret = new Accumulator(0.0001);
		int cnt = 0;
		for (double jeff = 0.0; jeff < 0.308; jeff += 0.004) {
			cnt++;
			IsingSolution sol = new IsingSolution(jeff);
			ret.accum(sol.J, sol.E0);
		}
		return ret;
	}
	
	public static void main(String[] args) {
		IsingSMF o = new IsingSMF();
		o.term = new Terminal();
		o.term.importObject(o);
		o.term.runApplication();
	}
}
