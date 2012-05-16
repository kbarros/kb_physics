package kip.javasim.quantum;

import kip.javasim.Random;
import scikit.numerics.fn.C1FunctionND;
import scikit.numerics.opt.Relaxation;
import scikit.util.DoubleArray;
import scikit.util.Pair;

public class IsingChain2 {
	int n = 6;
	int dim = n/2;
	double dt = 0.0001;
	
	Random r = new Random(1);
	
	void run() {
		Relaxation opt = new Relaxation(dim, dt);
		opt.setFunction(energyFn);

		double[] a = new double[dim]; // coeff
		for (int i = 0; i < dim; i++) {
			a[i] = r.nextDouble()-0.5;
//			a[i] = 0.1;
		}
		opt.initialize(a);
		opt.setStepSize(dt);
		
		while(true) {
			for (int i = 1; i < dim; i++)
				System.out.print(a[i] + " " );
			System.out.println();
			
			System.out.println(energyFn.eval(a));
			opt.step();
		}
	}
	
	C1FunctionND energyFn = new C1FunctionND() {
		double[] grad = new double[dim];
		
		double J = 1;
		double lambda = 100;
		
		double[] a;
		double[] spin = new double[n];
		double action; // a_ij s_i s_j
		
		double e_acc;
		double[] corr_acc = new double[dim];
		double[] e_corr_acc = new double[dim];
		
		
		public Pair<Double,double[]> calculate(final double[] a) {
			this.a = a;
			for (int i = 0; i < n; i++)
				spin[i] = -1;
			action = calcAction();
			
			e_acc = 0;
			DoubleArray.zero(corr_acc);
			DoubleArray.zero(e_corr_acc);
			
			boolean mc = false;
			int mcs = mc ? mcAccumulate() : enumAccumulate();
			
			double e = e_acc / mcs; 
			for (int i = 0; i < dim; i++) {
				double corr = corr_acc[i] / mcs;
				double e_corr = e_corr_acc[i] / mcs;
				grad[i] = 2*(e_corr - e*corr);
			}
			return new Pair<Double,double[]>(e, grad);
		}

		int mcAccumulate() {
			int mcs;
			for (mcs = 0; mcs < 1000; mcs++) {
				for (int i = 0; i < n; i++) {
					int j = r.nextInt(n);
					trySpinFlip(j);
				}
				accumulateConfig(1);
			}
			return mcs;
		}

		int enumAccumulate() {
			double norm = 0;
			for (int c = 0; c < 1<<n; c++) {
				for (int i = 0; i < n; i++)
					spin[i] = ((c & (1<<i)) == 0) ? 1 : -1;
				norm += Math.exp(2*calcAction());
			}
			
			for (int c = 0; c < 1<<n; c++) {
				for (int i = 0; i < n; i++)
					spin[i] = ((c & (1<<i)) == 0) ? 1 : -1;
				double weight = Math.exp(2*calcAction()) / norm;
				accumulateConfig(weight);
			}
			return 1;
		}
		
		double calcAction() {
			double action = 0;
			for (int i = 0; i < n; i++) {
				for (int d = 0; d < dim; d++) { 
					int ip = (i+d)%n;
					action += a[d] * spin[i] * spin[ip];
				}
			}
			return action;
		}
		
		double actionAfterFlip(int i) {
			double actionp = action;
			for (int j = i-dim+1; j < i+dim; j++) {
				if (i != j) {
					int d = Math.abs(j - i);
					int jp = (j+n)%n;
					actionp -= 2 * a[d] * spin[i] * spin[jp];
				}
			}
			
//			spin[i] *= -1;
//			assert(Math.abs(weightp - calcWeight()) < 1e-12);
//			spin[i] *= -1;
			
			return actionp;
		}
		
		void trySpinFlip(int i) {
			double actionp = actionAfterFlip(i);
			
			double p1 = Math.exp(2*action);
			double p2 = Math.exp(2*actionp);
			if (p2 > p1 || r.nextDouble() < p2/p1) {
				spin[i] *= -1;
				action = actionp;
			}
		}
		
		double classicalEnergy() {
			double energy = 0;
			for (int i = 0; i < n; i++) {
				int ip = (i+1)%n;
				energy += -J*spin[i]*spin[ip];
			}
			return energy;
		}
		
		double classicalCorrelation(int d) {
			double corr = 0;
			for (int i = 0; i < n; i++) {
				int ip = (i+d)%n;
				corr += spin[i]*spin[ip];
			}
			return corr/n;
		}
		
		void accumulateConfig(double w) {
			double e = classicalEnergy();
			for (int i = 0; i < n; i++) {
				double f = Math.exp(action);
				double fp = Math.exp(actionAfterFlip(i));
				e += lambda * fp / f;
			}
			
			e_acc += w*e;
			for (int d = 0; d < dim; d++) {
				double corr = classicalCorrelation(d);
				corr_acc[d] += w*corr;
				e_corr_acc[d] += w*e*corr;
			}
		}
	};

	
	public static void main(String[] args) {
		new IsingChain2().run();
	}
}
