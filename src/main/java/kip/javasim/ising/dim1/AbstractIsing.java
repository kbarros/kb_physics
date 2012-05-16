package kip.javasim.ising.dim1;

import scikit.jobs.params.Parameters;
import static java.lang.Math.*;


public abstract class AbstractIsing extends Dynamics1D {
	public enum DynType {METROPOLIS, GLAUBER, KAWA_GLAUBER, KAWA_METROPOLIS};
	public DynType dynamics = DynType.GLAUBER;
	public double T, J;

	public void initialize(Parameters params) {
		super.initialize(params);
		
		String dyn = params.sget("Dynamics", "Ising Glauber");
		if (dyn.equals("Ising Glauber"))
			dynamics = DynType.GLAUBER;
		else if (dyn.equals("Ising Metropolis"))
			dynamics = DynType.METROPOLIS;
		else if (dyn.equals("Kawasaki Glauber"))
			dynamics = DynType.KAWA_GLAUBER;
		else if (dyn.equals("Kawasaki Metropolis"))
			dynamics = DynType.KAWA_METROPOLIS;
		
		N = Integer.highestOneBit(params.iget("N"));
		params.set("N", N);
		
		R = min(params.iget("R"), N/2-1);
		params.set("R", R);
		
		dx = Integer.highestOneBit(params.iget("dx", 1));
		params.set("dx", dx);
	}
	
	public void setParameters(Parameters params) {
		super.setParameters(params);
		T  = params.fget("T");
		J  = params.fget("J", 1);
	}

}
