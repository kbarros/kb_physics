package kip.javasim.ising.dim2.apps;

import static java.lang.Math.exp;
import static java.lang.Math.min;
import static scikit.numerics.Math2.sqr;
import static scikit.util.Utilities.format;

import java.awt.Color;

import kip.javasim.Random;
import kip.javasim.ising.spinblock.SpinBlocks2D;
import scikit.dataset.Accumulator;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.Parameters;
import scikit.numerics.fft.FFT2D;


public class GeneralizedLinear2DApp extends Simulation {
	public static void main(String[] args) {
		new Control(new GeneralizedLinear2DApp(), "Ising Model");
	}
	
	Grid grid = new Grid("Coarse Grained Field");
	Grid structureGrid = new Grid("Structure Grid");
	Plot structurePlot = new Plot("Structure Evolution", "Time", "Structure");
	
	int dx;
	IsingLR sim;
	StructureFactor structure; 
	Accumulator structMax;
	
	public void load(Control c) {
		c.frame(grid);
		c.frame(structureGrid);
		c.frame(structurePlot);
		params.addm("Transition type", new ChoiceValue("Fluid->Clump", "Fluid->Stripe", "Stripe->Clump", "Clump->Fluid"));
		params.add("Random seed", 0);
		params.add("L", 512);
		params.add("R", 85);
		params.add("Initial magnetization", 0.0);
		params.addm("T", 0.05);
		params.addm("J", -1.0);
		params.addm("h", 0.8);
		params.addm("dt", 0.1);
		params.add("time");
		params.add("magnetization");
		
		flags.add("Res up");
		flags.add("Res down");
	}
	
	
	public void animate() {
		params.set("time", format(sim.t));
		params.set("magnetization", format(sim.magnetization()));
		sim.setParameters(params);
		
		int Lp = sim.L/dx;
		grid.registerData(Lp, Lp, sim.getField(dx));
		structureGrid.registerData(Lp, Lp, structure.data);
		structurePlot.registerPoints("Structure Max", structMax, Color.BLUE);
		
		if (flags.contains("Res up"))
			sim.scaleUp(params);
		if (flags.contains("Res down"))
			sim.scaleDown(params);
		flags.clear();
	}
	
	public void clear() {
		grid.clear();
		structureGrid.clear();
	}
	
	public void run() {
		grid.setScale(-1, 1);
		
		sim = new IsingLR(params);
		sim.randomizeField(params.fget("Initial magnetization"));
		dx = Math.max(Integer.highestOneBit(sim.R)/8, 1);
		structure = new StructureFactor();
		structMax = new Accumulator();
		
		while (true) {
			equilibrate();
			quenchParams();
			sim.t = 0;
			
			int dominantK = getGrowthMode();
			System.out.println("dom " + dominantK);
			
			while (sim.t < 10) {
				structure.calculate();
				structMax.accum(sim.t, structure.data[dominantK]);
				sim.step();
				
				Job.animate();
			}
		}
	}
	
	void equilibrate() {
		String type = params.sget("Transition type");
		if (type.equals("Fluid->Clump") || type.equals("Fluid->Stripe")) {
			int L = sim.L;
			for (int i = 0; i < L*L; i++) {
				sim.spins.set(i % L, i / L, sim.random.nextDouble() < 0.5 ? -1 : 1);
			}
			sim.t = 0;
		}
		else {
			if (type.equals("Stripe->Clump")) {
				params.set("h", 0.0);
				params.set("T", 0.05);
			}
			if (type.equals("Clump->Fluid")) {
				params.set("h", 0.8);
				params.set("T", 0.05);			
			}
			sim.setParameters(params);
			sim.t = 0;
			while (sim.t < 20) {
				sim.step();
				// structure.calculate();
				Job.animate();
			}
			sim.t = 0;
		}
	}
	
	void quenchParams() {
		String type = params.sget("Transition type");
		if (type.equals("Fluid->Clump") || type.equals("Stripe->Clump")) {
			params.set("h", 0.8);
			params.set("T", 0.05);			
		}
		if (type.equals("Fluid->Stripe")) {
			params.set("h", 0.0);
			params.set("T", 0.05);
			
		}
		if (type.equals("Clump->Fluid")) {
			params.set("h", 0);
			params.set("T", 5);
		}
		sim.setParameters(params);
	}
	
	int getGrowthMode() {
		int k0 = 2; // valid when L=256, R=85
		String type = params.sget("Transition type");
		if (type.equals("Fluid->Clump") || type.equals("Fluid->Stripe") || type.equals("Clump->Fluid"))
			return k0;
		if (type.equals("Stripe->Clump")) {
			structure.calculate();
			double[] d = structure.data;
			int Lp = sim.L/dx;
			return (d[k0] < d[k0*Lp]) ? k0 : k0*Lp;
		}
		return -1;
	}
	
	class StructureFactor {
		FFT2D fft;
		double[] data;
		
		public StructureFactor() {
			int Lp = sim.L / dx;
			fft = new FFT2D((int)(Lp), (int)(Lp));
			fft.setLengths(sim.L, sim.L);
			
			data = new double[Lp*Lp];
		}
		
		public void calculate() {
			int Lp = fft.dim1;
			
			double[] scratch = fft.getScratch();
			fft.transform(sim.getField(dx), scratch);
			
			double vol = sim.L*sim.L;
			for (int i = 0; i < Lp*Lp; i++) {
				data[i] = (sqr(scratch[2*i+0])+sqr(scratch[2*i+1]))/vol; 
			}
		}
	}

	
	class IsingLR {
		public SpinBlocks2D spins;
		public double dt, t;
		public int L, R;
		public double T, J, h;
		public Random random = new Random();
		
		public IsingLR(Parameters params) {
			L = Integer.highestOneBit(params.iget("L"));
			params.set("L", L);
			R = min(params.iget("R"), L/2-1);
			params.set("R", R);
			
			spins = new SpinBlocks2D(L, R);
			random.setSeed(params.iget("Random seed"));
			setParameters(params);
			
			t = 0;
		}
		
		public void setParameters(Parameters params) {
			dt = params.fget("dt");
			T  = params.fget("T");
			J  = params.fget("J", 1);
			h  = params.fget("h", 0);
		}		
		
		public void scaleUp(Parameters params) {
			SpinBlocks2D newSpins = new SpinBlocks2D(2*L, 2*R);
			for (int x = 0; x < L; x++) {
				for (int y = 0; y < L; y++) {
					int s = spins.get(x, y);
					newSpins.set(2*x+0, 2*y+0, s);
					newSpins.set(2*x+0, 2*y+1, s);
					newSpins.set(2*x+1, 2*y+0, s);
					newSpins.set(2*x+1, 2*y+1, s);
				}
			}
			spins = newSpins;
			L *= 2;
			R *= 2;
			params.set("L", L);
			params.set("R", R);
		}
		
		
		public void scaleDown(Parameters params) {
			L /= 2;
			R /= 2;
			SpinBlocks2D newSpins = new SpinBlocks2D(L, R);
			for (int x = 0; x < L; x++) {
				for (int y = 0; y < L; y++) {
					int s = spins.get(2*x, 2*y);
					newSpins.set(x, y, s);
				}
			}
			spins = newSpins;
			params.set("L", L);
			params.set("R", R);
		}
		

		public double magnetization() {
			return (double)spins.sumAll() / (L*L);
		}

		public void randomizeField(double m) {
			for (int i = 0; i < L*L; i++) {
				// p(s = +-1) = (1 +- m) / 2
				int s = (random.nextDouble() < (1+m)/2) ? 1 : -1;
				spins.set(i%L, i/L, s);
			}
		}
		
		public double[] getField(int dx) {
			int scale = Integer.numberOfTrailingZeros(dx);
			int blocks[] = spins.blocksAtScale(scale);
			double ret[] = new double[blocks.length];
			double blockSize = (1<<scale)*(1<<scale);
			for (int i = 0; i < blocks.length; i++)
				ret[i] = blocks[i]/blockSize;
			return ret;
		}
		
		
		private boolean shouldFlip(double dE) {
//			return dE <= 0 || random.nextDouble() < Math.exp(-dE/T);
			return (dE <= 0 && T == 0) || random.nextDouble() < exp(-dE/T)/(1+exp(-dE/T));
		}
		
		protected void step() {
			for (int cnt = 0; cnt < L*L*dt; cnt++) {
				int x1 = random.nextInt(L);
				int y1 = random.nextInt(L);
				int s1 = spins.get(x1, y1);
				double dE = 2*s1*(h + J*(spins.sumInRange(x1,y1)-s1)/(4*R*R));
				if (shouldFlip(dE)) {
					spins.flip(x1, y1);
				}
				scikit.jobs.Job.yield();
			}
			t += dt;
		}
	}
}
