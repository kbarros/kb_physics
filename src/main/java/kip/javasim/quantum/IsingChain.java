package kip.javasim.quantum;

public class IsingChain {
	int n = 5;
	int mask = (1 << n) - 1;
	double[] ground;
	double[] ground_deriv;
	double g = 1; // coupling
	
	public IsingChain() {
		ground = new double[1<<n];
		ground_deriv = new double[1<<n];
		
//		for (int c = 0; c < 1<<n; c++)
//			ground[c] = Math.random()-0.5;
		
		for (int c = 0; c < 1<<n; c++)
			ground[c] = 1;
	}
	
	double sqr(double x) {
		return x*x;
	}
	
	int rotateOne(int c) {
		c = c & mask;
		return ((c << 1) & mask) | (c >> (n-1));
	}
	
	double classicalEnergy(int c) {
		int cnt = Integer.bitCount(c ^ rotateOne(c));
		return 2*cnt - n;
	}
	
	double hamiltonian(int c1, int c2) {
		if (c1 == c2)
			return classicalEnergy(c1);
		else if (Integer.bitCount(c1 ^ c2) == 1) 
			return g;
		else
			return 0;
	}
	
	double trialHamiltonian(int c) {
		double ret = 0;
		for (int i = 0; i < n; i++) {
			ret += ground[c ^ (1 << i)];
		}
		return -g * ret / ground[c];
	}
	
	double energyEstimate() {
		double dh_mean = 0;
		for (int c = 0; c < (1<<n); c++) {
			double hp = trialHamiltonian(c);
			double h = hamiltonian(c,c);
			dh_mean += h - hp;
		}
		return dh_mean / (1<<n);
	}
	
	double hamiltonianDistance() {
		double dh_mean = energyEstimate();
		double ret = 0;
		for (int c = 0; c < (1<<n); c++) {
			double hp = trialHamiltonian(c);
			double h = hamiltonian(c,c);
			ret += sqr((h - hp) - dh_mean);
		}
		return ret / (1<<n);
	}
	
	void calc_grad() {
		double dh_mean = energyEstimate();
		
		for (int c = 0; c < (1<<n); c++) {
			// d <(Del H)^2> / d A_c
			ground_deriv[c] = 0;
			double gs_sum = 0;
			
			// contribution from neighbors
			for (int i = 0; i < n; i++) {
				int c2 = (c ^ (1 << i));
				double hp = trialHamiltonian(c2);
				double h = hamiltonian(c2,c2);
				ground_deriv[c] += -2*g*(h - hp - dh_mean) / ground[c2];
				gs_sum += ground[c2];
			}
			// contribution from self term
			double hp = trialHamiltonian(c);
			double h = hamiltonian(c,c);
			ground_deriv[c] += 2*g*(h - hp - dh_mean) * gs_sum / sqr(ground[c]);
		}
	}
	
	void step() {
		calc_grad();
		
		double dt = 0.0001;
		for (int c = 0; c < (1<<n); c++) {
			ground[c] += dt * ground_deriv[c];
		}
	}
	
	void printHamiltonian() {
		for (int c1 = 0; c1 < (1<<n); c1++) {
			for (int c2 = 0; c2 < (1<<n); c2++) {
				System.out.print(hamiltonian(c1, c2) + "    ");
			}
			System.out.println();
		}
	}
	
	void printHamiltonianDistance() { 
		System.out.println(energyEstimate());
		System.out.println(hamiltonianDistance());
	}
	
	public static void main(String[] args) {
		IsingChain chain = new IsingChain();
		
		for (int i = 0; i < 10000; i++) {
			chain.step();
		}
		
		chain.printHamiltonianDistance();
	}
}
