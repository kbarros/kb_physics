package kip.javasim.clump.dim2;

import static scikit.numerics.Math2.j0;
import static scikit.numerics.Math2.j1;
import static scikit.numerics.Math2.jn;
import kip.javasim.Random;
import scikit.dataset.Accumulator;
import scikit.jobs.params.Parameters;

abstract public class AbstractClump2D {
	public double L, R, T, dx;
	Random random = new Random();
	
	public static final double DENSITY = 1;
	
	// value of kR which minimizes j1(kR)/kR
	public static final double KR_SP = 5.13562230184068255630140;
	// S(k) ~ 1 / (V(kR_sp)/T+1)
	// => T_SP = - V(kR_sp) = - 2 j1(kR_sp) / kR_sp 
	public static final double T_SP = 0.132279487396100031736846;	
	
	public double potential(double kR) {
		return (kR == 0 ? 1 : 2*j1(kR)/kR);
	}
	
	public double dpotential_dkR(double kR) {
		double kR2 = kR*kR;
		return (kR == 0) ? 0 : j0(kR)/kR - 2*j1(kR)/kR2  - jn(2,kR)/kR;
	}
	
	abstract public void readParams(Parameters params);
	abstract public Accumulator newStructureAccumulator(double binWidth);
	abstract public void accumulateStructure(Accumulator sf);
	abstract public void simulate();
	abstract public double[] coarseGrained();
	abstract public int numColumns();
	abstract public double time();
}
