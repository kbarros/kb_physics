package kip.javasim;

public class AutoCorrelation {
	int[] a0;
	double mean, var;
	
	public AutoCorrelation(int[] data) {
		a0 = new int[data.length];
		System.arraycopy(data, 0, a0, 0, data.length);
		
		mean = var = 0;
		for (int i = 0; i < a0.length; i++) {
			mean += a0[i];
			var += a0[i]*a0[i];
		}
		mean /= a0.length;
		var /= a0.length;
	}
	
	public double calc(int[] a1) {
		double acc = 0;
		for (int i = 0; i < a0.length; i++) {
			acc += (a0[i]-mean)*(a1[i]-mean);
		}
		return acc / var;
	}
}
