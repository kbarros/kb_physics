package kip.javasim.fun;


import static java.lang.Math.PI;
import static java.lang.Math.exp;
import static java.lang.Math.random;
import static java.lang.Math.sqrt;

import java.awt.Color;

import scikit.dataset.Accumulator;
import scikit.dataset.Function;
import scikit.dataset.Histogram;
import scikit.dataset.PointSet;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DoubleValue;
import scikit.jobs.params.IntValue;


public class PistonApp extends Simulation {
	public static void main(String[] args) {
		new Control(new PistonApp(), "Ideal Gas Simulation");
	}
	
	double[] pistonLine = new double[2];
	Plot particles = new Plot("Particles");
	Plot enthalpy = new Plot("Enthalpy: U + PV");
	Plot idealGas = new Plot("PV / NkT");
	Plot distrib = new Plot("Velocity distribution"); // 0.01
	Accumulator idealGasAcc, kineticAcc, workAcc, energyAcc;
	
	int N;			// number of particles
	double T;		// kinetic energy
	double t, dt;	// time, time step
	
	double px, pv;	// piston phase coordinates
	double pa;		// piston acceleration
	double pm;		// piston mass
	double[] x, v;	// gas phase coordinates
	double m = 1;	// gas mass
	
	
	public void load(Control c) {
		c.frameTogether("Plots", particles, enthalpy, idealGas, distrib);
		params.add("Initial piston position", new DoubleValue(10.0, 0, 100));
		params.add("Initial piston velocity", new DoubleValue(0.0, -1, 1));
		params.add("# of particles", new IntValue(1000, 1, 5000));
		params.addm("Piston mass", new DoubleValue(100.0, 0, 200).withSlider());
		params.addm("Piston acceleration", new DoubleValue(0.00001, 0,  0.00002).withSlider());
		params.addm("dt", new DoubleValue(0.05, 0, 0.2));
		params.addm("Bin width", new DoubleValue(0.0002, 0.00005, 0.01));
	}
	
	
	public void animate() {
		pm  = params.fget("Piston mass");
		pa	= params.fget("Piston acceleration");
		dt	= params.fget("dt");

		pistonLine[0] = pistonLine[1] = px;
		particles.registerPoints("Points", new PointSet(0, 1, x), Color.BLACK);
		particles.registerLines("", new PointSet(-0.1, N-0.8, pistonLine), Color.BLUE);
		
		idealGas.registerLines("Ideal gas", idealGasAcc, Color.BLACK);
		
		enthalpy.registerLines("Kinetic", kineticAcc, Color.BLACK);
		enthalpy.registerLines("Work", workAcc, Color.BLUE);
		enthalpy.registerLines("Energy", energyAcc, Color.RED);
		
		Histogram velocities = new Histogram(params.fget("Bin width"));
		velocities.setNormalizing(true);
		for (double vi : v)
			velocities.accum(vi);
		distrib.registerBars("Velocity distribution", velocities, Color.RED);
		distrib.registerLines("Gaussian distribution", new Function() {
			final double beta = 1 / kT();
			public double eval(double v) {
				return velocityProbability(v, beta);
			}
		}, Color.BLACK);
	}
	
	public void clear() {
		particles.clear();
		idealGas.clear();
		enthalpy.clear();
		distrib.clear();
	}
	
	public void run() {		
		N = params.iget("# of particles");
		px	= params.fget("Initial piston position");
		pv	= params.fget("Initial piston velocity");
		pm  = params.fget("Piston mass");
		pa	= params.fget("Piston acceleration");
		dt	= params.fget("dt");
		
		double bw = 5000*dt;
		idealGasAcc = new Accumulator(bw);
		kineticAcc = new Accumulator(bw);
		workAcc = new Accumulator(bw);
		energyAcc = new Accumulator(bw);
		
		t = 0;
		v = new double[N];
		x = new double[N];
		for (int i = 0; i < N; i++) {
			x[i] = random() * px;
			v[i] = 0.001 * (2*random() - 1);
		}
		
		while (true) {
			simulationStep();
			
			idealGasAcc.accum(t, force()*px / (N * kT()));
			kineticAcc.accum(t, kineticEnergy());
			workAcc.accum(t, force()*px);
			energyAcc.accum(t, kineticEnergy() + force()*px);

			Job.animate();
		}
	}
	
	
	private void simulationStep() {
		pv += - pa * dt;
		px += pv * dt - 0.5*pa*dt*dt;
		
		for (int i = 0; i < N; i++) {
			x[i] += v[i] * dt;
			
			if (x[i] < 0 && v[i] < 0) {
				x[i] = 0;
				v[i] = -v[i];
			}
			
			if (x[i] > px && v[i] > pv) {
				// elastic collision with piston
				x[i] = px;
				double m1 = m;
				double m2 = pm;
				double v1 = v[i];
				double v2 = pv;
				v[i] =   (m1-m2)*v1/(m1+m2) + 2*m2*v2/(m1+m2);
				pv   = - (m1-m2)*v2/(m1+m2) + 2*m1*v1/(m1+m2);
			}
		}
		
		t += dt;
	}
	
	// F = ma
	private double force() {
		return pa*pm;
	}
	
	// kinetic energy = internal energy = 1/2 NkT
	private double kT() {
		return 2*kineticEnergy()/N;
	}
	
	private double kineticEnergy() {
		double sum = 0;
		for (double vi : v)
			sum += vi*vi;
		return sum/2;
	}
	
	// Probability density P particle velocity lies within [v, v+dv]
	// Boltzmann velocity distribution: P*dv ~ e^(beta*m*v^2 / 2)
	// Normalize by integral[P*dv, v=-inf..inf] = sqrt(2*pi / beta*m)
	private double velocityProbability(double v, double beta) {
		return sqrt(beta*m/(2*PI))*exp(-beta*m*v*v/2);
	}
}
