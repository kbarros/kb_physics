package kip.javasim.ising.dim1.apps;

import static java.lang.Math.PI;
import static java.lang.Math.exp;
import static java.lang.Math.sin;
import static scikit.util.Utilities.format;

import java.awt.Color;

import kip.javasim.ising.dim1.AbstractIsing;
import kip.javasim.ising.dim1.FieldIsing;
import kip.javasim.ising.dim1.Ising;
import scikit.dataset.Accumulator;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;


// TODO redundant
class Structure {
	int Lp;
	double L;
	double R;
	double kRmin = 0, kRmax = Double.MAX_VALUE;
	
	scikit.numerics.fft.managed.RealDoubleFFT fft;
	private double[] fftScratch;
	
	public Structure(int Lp, double L, int R) {
		this.Lp = Lp;
		this.L = L;
		this.R = R;
		fftScratch = new double[Lp];
		fft = new scikit.numerics.fft.managed.RealDoubleFFT_Radix2(Lp);
	}

	public void setBounds(double kRmin, double kRmax) {
		this.kRmin = kRmin;
		this.kRmax = kRmax;
	}
	
	public void accumulate(double[] field, Accumulator acc) {
		double dx = L/Lp;
		for (int i = 0; i < Lp; i++)
			fftScratch[i] = field[i];
		fft.transform(fftScratch);
		
		for (int i = 0; i < Lp/2; i++) {
			double kR = 2*PI*i*R/L;
			if (kR >= kRmin && kR <= kRmax) {			
				double re = fftScratch[i];
				double im = (i == 0) ? 0 : fftScratch[Lp-i];
				acc.accum(kR, (re*re+im*im)*dx*dx/L);
			}
		}
	}
	
	public double theory(AbstractIsing sim, double kR) {
		double Q = kR == 0 ? 1 : sin(kR)/kR;
		double K = sim.J/sim.T;
		double t = sim.time();
		double M = 2;
		double D;
		
		switch (sim.dynamics) {
		case METROPOLIS:
			M = 4; // fall through
		case GLAUBER:
			D = -1 + Q*K;
			return (exp(M*D*t)*(1 + 1/D) - 1/D);
			
		case KAWA_METROPOLIS:
			M = 4; // fall through
		case KAWA_GLAUBER:
			D = -1 + Q*(1 + K) - Q*Q*K;
			return  kR == 0 ? 1 : exp(M*D*t)*(1+(1-Q)/D) - (1-Q)/D;

		default:
			return Double.NaN;
		}
	}
	
	
	// TODO use function
	public Accumulator coarseGrainedTheory(AbstractIsing sim, double binWidth) {
		Accumulator ret = new Accumulator(binWidth);
		for (int i = 0; i < Lp/2; i++) {
			double kR = 2*PI*i*R/L;
			if (kR >= kRmin && kR <= kRmax)
				ret.accum(kR, theory(sim, kR));
		}
		return ret;
	}
}


public class OrderingApp extends Simulation {
	Plot structurePlot = new Plot("Structure");
	AbstractIsing sim;
	Structure structure;
	double binWidth;
	Accumulator[] structData;
	double[] field;
	int stepCnt; 
	int numSteps = 10;
	
	public static void main(String[] args) {
		new Control(new OrderingApp(), "Growth for Ising Droplets");
	}

	public void load(Control c) {
		c.frame(structurePlot);
		params.add("Dynamics", new ChoiceValue("Ising Glauber", "Ising Metropolis", "Kawasaki Glauber", "Kawasaki Metropolis"));
		params.add("Simulation type", new ChoiceValue("Ising", "Langevin"));
		params.add("kR maximum", 20.0);
		params.add("kR bin width", 0.1);
		params.add("Random seed", 0);
		params.add("N", 1<<20);
		params.add("R", 1<<12);
		params.add("dx", 1<<6);
		params.add("T", 4.0/9.0);
		params.add("J", 1.0);
		params.add("dt", 0.1);
		params.add("time");
	}
	
	public void animate() {
		params.set("time", format(sim.time()));
		Accumulator structTheory = structure.coarseGrainedTheory(sim, binWidth);
		structurePlot.registerLines("Structure data", structData[stepCnt], Color.BLACK);
		structurePlot.registerLines("Structure theory", structTheory, Color.BLUE);
//		structurePlot.registerLines("Comparison", new Function() {
//			public double eval(double kR) {
//				return (structTheory.eval(kR)-structData[stepCnt].eval(kR))*sqrt(sim.R);
//			}
//		}, Color.RED);
	}
	
	public void clear() {
		structurePlot.clear();
	}
	
	public void run() {
		String type = params.sget("Simulation type");
		sim = type.equals("Ising") ? new Ising(params) : new FieldIsing(params);
		structure = new Structure(sim.N/sim.dx, sim.N, sim.R);
		structure.setBounds(0, params.fget("kR maximum"));
		structData = new Accumulator[numSteps];
		binWidth = params.fget("kR bin width");
		for (int i = 0; i < numSteps; i++) {
			structData[i] = new Accumulator(binWidth);
		}
		
		while (true) {
			sim.initialize(params);
			sim.randomizeField(0);
//			sim.setField(0);
			
			for (stepCnt = 0; stepCnt < numSteps; stepCnt++) {
				sim.step();
				structure.accumulate(sim.copyField(), structData[stepCnt]);				
				Job.animate();
			}
			
			params.set("Random seed", params.iget("Random seed")+1);			
		}
	}
}
