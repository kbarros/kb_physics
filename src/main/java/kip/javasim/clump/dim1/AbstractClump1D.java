package kip.javasim.clump.dim1;

import static java.lang.Math.cos;
import static java.lang.Math.sin;
import kip.javasim.Random;
import scikit.dataset.Accumulator;
import scikit.jobs.params.Parameters;

abstract public class AbstractClump1D {
	public double L, R, T, dx;
	Random random = new Random();

	public static final double DENSITY = 1;
	public static final double KR_SP = 4.4934094579090641753;
	public static final double T_SP = 0.21723362821122165741;
	
	
	public double potential(double kR) {
		return (kR == 0) ? 1 : sin(kR)/kR;
	}
	
	public double dpotential_dkR(double kR) {
		return (kR == 0) ? 0 : -sin(kR)/(kR*kR) + cos(kR)/kR;
	}
	
	abstract public void readParams(Parameters params);
	abstract public Accumulator newStructureAccumulator(double binWidth);
	abstract public void accumulateStructure(Accumulator sf);
	abstract public void simulate();
	abstract public double[] coarseGrained();
	abstract public int numColumns();
	abstract public double time();
}
