package kip.javasim.ising.dim1;

import scikit.jobs.Job;
import scikit.jobs.params.Parameters;
import kip.javasim.ising.spinblock.SpinBlocks1D;
import static java.lang.Math.*;


public class Ising extends AbstractIsing {
	public SpinBlocks1D spins;
	
	
	public Ising(Parameters params) {
		initialize(params);
	}
	
	
    public Ising clone() {
		Ising c = (Ising)super.clone();
		c.spins = (SpinBlocks1D)spins.clone();
		return c;
    }
	
	
	public void initialize(Parameters params) {
		super.initialize(params);
		spins = new SpinBlocks1D(N, R, -1);
	}
	
	
	public double magnetization() {
		return (double)spins.sumAll() / N;
	}
	

	public void randomizeField(double m) {
		if (m == 1 || m == -1) {
			for (int i = 0; i < N; i++)
				spins.set(i, (int)m);
		}
		else {
			for (int i = 0; i < N; i++) {
				// p(s = +-1) = (1 +- m) / 2
				int s = (random.nextDouble() < (1+m)/2) ? 1 : -1;
				spins.set(i, s);
				Job.yield();
			}
		}
	}
	
	
	public void setField(double m) {
		double mAcc = 0;
		for (int i = 0; i < N; i++) {
			int s = (mAcc > i*m) ? -1 : 1;
			spins.set(i, s);
			mAcc += s;
		}
	}
	
	public double fieldElement(int i) {
		int scale = Integer.numberOfTrailingZeros(dx);
		return spins.getBlock(scale,i)/(double)dx;
	}
	
	
	private boolean shouldFlip(double dE) {
		switch (dynamics) {
			case METROPOLIS:
			case KAWA_METROPOLIS:
				return dE <= 0 || random.nextDouble() < Math.exp(-dE/T);
			case GLAUBER:
			case KAWA_GLAUBER:
				return random.nextDouble() < exp(-dE/T)/(1+exp(-dE/T));
			default:
				assert false;
		}
		return false;
	}
	
	protected void _step() {
		for (int cnt = 0; cnt < N*dt; cnt++) {
			int i = random.nextInt(N);
			int s_i = spins.get(i);
			double dE = 2*s_i*(h + J*(spins.sumInRange(i)-s_i)/(2*R));
			scikit.jobs.Job.yield();
			switch (dynamics) {
				case METROPOLIS:
				case GLAUBER:
					if (shouldFlip(dE))
						spins.flip(i);
					break;
					
				case KAWA_GLAUBER:
				case KAWA_METROPOLIS:
					int j = i + random.nextInt(2*R+1)-R;
					j = (j + N) % N;
					int s_j = spins.get(j);
					if (s_j != s_i) {
						dE += 2*s_j*(h + J*(spins.sumInRange(j)-s_j)/(2*R));
						if (shouldFlip(dE)) {
							spins.flip(i);
							spins.flip(j);
						}
					}
					break;
				default:
					assert false;
			}
		}
	}
}
