package kip.javasim.langevin.apps;

import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.util.DoubleArray;
import scikit.util.Utilities;
import static scikit.numerics.Math2.sqr;

public class Phi4App extends Simulation {
	Grid grid = new Grid("Grid");
	Phi4 sim;
	
	public static void main(String args[]) {
		new Control(new Phi4App(), "Phi4 Dynamics");
	}
	
	public void animate() {
		grid.registerData(sim.Lp, sim.Lp, sim.phi);
		sim.R = params.fget("R");
		sim.T = params.fget("T");
		sim.m = params.fget("m");
		sim.dt = params.fget("dt");
		
		params.set("free energy", Utilities.format(sim.freeEnergy() / sqr(sim.L)));
		
		if (flags.contains("Bomb"))
			sim.bomb();
		if (flags.contains("Clean"))
			sim.clean();
		flags.clear();
	}

	public void clear() {
		grid.clear();
	}

	public void load(Control c) {
		c.frame(grid);
		params.add("L", 100.0);
		params.add("dx", 1.0);
		params.addm("R", 1.0);
		params.addm("T", -1.0);
		params.addm("m", -0.8);
		params.addm("dt", 0.05);
		
		params.add("free energy");
		
		flags.add("Bomb");
		flags.add("Clean");
	}

	public void run() {
		sim = new Phi4(params.fget("L"), params.fget("dx"));
		sim.R = params.fget("R");
		sim.T = params.fget("T");
		sim.m = params.fget("m");
		sim.dt = params.fget("dt");
		
		while (true) {
			sim.step();
			Job.animate();
		}
	}
	
}

class Phi4 {
	double L, dx;
	double R, T, m, dt;
	int Lp;
	double[] phi;
	double[] scratch;
	
	public Phi4(double L, double dx) {
		this.L = L;
		this.dx = dx;
		Lp = (int)(L/dx);
		phi = new double[Lp*Lp];
		scratch = new double[Lp*Lp];
		
		for (int i = 0; i < Lp*Lp; i++) {
			phi[i] = 0.01*(Math.random() - 0.5);
		}
	}
	
	// F = \int dx [ (1/2) (R \del phi)^2 + (1/4) (phi^2 + T)^2 ]  
	public double freeEnergy() {
		double ret = 0;
		for (int i = 0; i < Lp*Lp; i++) {
			int il = (i % Lp == 0)    ? i+(Lp-1) : i-1;
			int ir = (i % Lp == Lp-1) ? i-(Lp-1) : i+1;
			int id = (i / Lp == 0)    ? i+Lp*(Lp-1) : i-Lp;
			int iu = (i / Lp == Lp-1) ? i-Lp*(Lp-1) : i+Lp;
			
			// (\del phi)^2 term
			ret += 0.5 * sqr(R*(phi[ir]-phi[il])/(2*dx));
			ret += 0.5 * sqr(R*(phi[iu]-phi[id])/(2*dx));
			
			// driving term
			ret += 0.25 * sqr(sqr(phi[i]) + T);
		}
		return ret*sqr(dx);
	}
	
	public void step() {
		for (int i = 0; i < Lp*Lp; i++) {
			int il = (i % Lp == 0)    ? i+(Lp-1) : i-1;
			int ir = (i % Lp == Lp-1) ? i-(Lp-1) : i+1;
			int id = (i / Lp == 0)    ? i+Lp*(Lp-1) : i-Lp;
			int iu = (i / Lp == Lp-1) ? i-Lp*(Lp-1) : i+Lp;
			
			scratch[i] = 0;
			
			// \del^2 phi term
			scratch[i] += - sqr(R)*(phi[ir]-2*phi[i]+phi[il])/sqr(dx);
			scratch[i] += - sqr(R)*(phi[iu]-2*phi[i]+phi[id])/sqr(dx);
			
			// driving term
			scratch[i] += (sqr(phi[i]) + T) * phi[i];
		}
		
		for (int i = 0; i < Lp*Lp; i++) {
			phi[i] += - dt * scratch[i];
		}
		
		// phi += m - <phi>
		DoubleArray.shift(phi, m-DoubleArray.mean(phi));
	}
	
	public void bomb() {
		for (int i = 0; i < Lp*Lp; i++)
			phi[i] = -Math.sqrt(-T);
		
		for (int i = 0; i < Lp/2; i++) {
			for (int j = 0; j < Lp/2; j++) {
				phi[i*Lp+j] = +Math.sqrt(-T); 
			}
		}
	}
	
	public void clean() {
		for (int i = 0; i < Lp*Lp; i++)
			phi[i] = 0;
	}

}
