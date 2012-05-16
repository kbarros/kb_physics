package kip.javasim.clump.dim2.apps;

import static java.lang.Math.PI;
import static java.lang.Math.cos;
import static java.lang.Math.log;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;
import static scikit.numerics.Math2.hypot;
import static scikit.numerics.Math2.j1;
import static scikit.numerics.Math2.sqr;
import static scikit.util.Utilities.format;
import kip.javasim.Random;
import scikit.graphics.GrayScale;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.numerics.fft.FFT2D;
import scikit.numerics.fn.C1FunctionND;
import scikit.numerics.fn.Function2D;
import scikit.numerics.opt.Constraint;
import scikit.numerics.opt.Relaxation;
import scikit.util.DoubleArray;
import scikit.util.Pair;

public class Saddle2App extends Simulation {
	final double inf = Double.POSITIVE_INFINITY;
	
	Grid grid = new Grid("Grid");
	Plot plot = new Plot("");
	
	double L, R, T, dt;
	int dim;
	double[] phi, phibar;
	Random random = new Random();
	FFT2D fft;
	double fe;
	int time;
	double targetVariance;
	
	
	public static void main(String[] args) {
		new Control(new Saddle2App(), "Clump Optimizer");
	}

	public void load(Control c) {
		c.frame(grid, plot);
		params.addm("Zoom", new ChoiceValue("Yes", "No"));
		params.addm("T", 0.14);
		params.addm("dt", 0.1);
		params.addm("var", 0.01);
		params.add("R", 1000.0);
		params.add("L", 20000.0);
		params.add("dim", 64);
		params.add("Random seed", 0);
		params.add("Time");
		params.add("F density");
	}
	
	public void animate() {
		T = params.fget("T");
		dt = params.fget("dt");
		
		grid.setColors(new GrayScale());
		if (params.sget("Zoom").equals("Yes"))
			grid.setAutoScale();
		else
			grid.setScale(0, 4);
		grid.registerData(dim, dim, phi);
		
		params.set("var", format(targetVariance)); 
//		targetVariance = params.fget("var");
		
		params.set("F density", format(fe-0.5));
		params.set("Time", time);
//		plot.registerLines("", new PointSet(0, 1, section), Color.BLUE);
	}
	
	public void clear() {
		grid.clear();
		plot.clear();
	}
	
	public void run() {
		T = params.fget("T");
		dt = params.fget("dt");
		R = params.fget("R");
		L = params.fget("L");
		dim = params.iget("dim");
		targetVariance = params.fget("var");
		time = 0;
		
		phi = new double[dim*dim];
		phibar = new double[dim*dim];
		initializeField();
		
		fft = new FFT2D(dim, dim);
		fft.setLengths(L, L);
		
		random.setSeed(params.iget("Random seed"));
		
		C1FunctionND freeEnergy = new C1FunctionND() {
			double[] grad = new double[dim*dim];
			public Pair<Double,double[]> calculate(final double[] p) {
				return calculateFreeEnergy(p, grad);
			}
		};
		
		Constraint varianceCst = new Constraint() {
			double[] grad = new double[dim*dim];
			public double stiffness() { return 0.2; }
			public Pair<Double,double[]> calculate(double[] p) {
				return calculateVariance(p, grad);
			}
		};
		
		Relaxation opt = new Relaxation(dim*dim, dt);
		opt.setFunction(freeEnergy);
		opt.addConstraint(varianceCst);
		opt.initialize(phi);
		
		while(true) {
			opt.setStepSize(dt);
			opt.step();
			targetVariance -= dt * opt.dc_dt(phi, varianceCst);
			fe = freeEnergy.eval(phi);
			System.out.println("++ " + opt.dc_dt(phi, varianceCst) + " " + varianceCst.eval(phi));
			time++;
			Job.animate();
		}
	}
	
	private Function2D potential = new Function2D() {
		public double eval(double k1, double k2) {
			double kR = hypot(k1,k2)*R;
			return kR == 0 ? 1 : 2*j1(kR)/kR;
		}
	};
	
	private void initializeField() {
		for (int i = 0; i < dim*dim; i++) {
			double dx = L/dim;
			double x = dx*(i%dim - dim/2);
			double y = dx*(i/dim - dim/2);
			double r = sqrt(x*x+y*y);
			double mag = 0.8 / (1+sqr(r/R));
			
			double KR_SP = 5.13562230184068255630140;
			double kR = KR_SP; // it's fun to try different values
			double x1 = x*cos(1*PI/6) + y*sin(1*PI/6);
			double x2 = x*cos(3*PI/6) + y*sin(3*PI/6);
			double x3 = x*cos(5*PI/6) + y*sin(5*PI/6);
			phi[i] = 1+mag*(cos(x1*kR/R) + cos(x2*kR/R) + cos(x3*kR/R));
//			phi[i] = 1+mag*random.nextGaussian()/5;	// random initial condition
		}
		
		double mean = DoubleArray.mean(phi);
		for (int i = 0; i < dim*dim; i++)
			phi[i] += 1-mean;
	}
	
	public Pair<Double,double[]> calculateFreeEnergy(double[] p, double[] grad) {
		final double[] pb = phibar;
		fft.convolve(p, pb, potential);
		double fe_acc = 0;
		for (int i = 0; i < p.length; i++) {
			if (p[i] <= 0) {
				fe_acc = grad[i] = inf;
			}
			else {
				fe_acc += (p[i]*pb[i]/2+T*p[i]*log(p[i])) / p.length;
				grad[i] = pb[i]+T*log(p[i]);
			}
		}
		double mu = DoubleArray.mean(grad);
		for (int i = 0; i < p.length; i++) {
			grad[i] -= mu;
		}
		return new Pair<Double,double[]>(fe_acc, grad);
	}

	public Pair<Double,double[]> calculateVariance(double[] p, double[] grad) {
		double c = 0;
		for (int i = 0; i < p.length; i++) {
			c += (p[i]-1)*(p[i]-1) / p.length;
			grad[i] = 2*(p[i]-1);
		}
		c -= targetVariance;
		return new Pair<Double,double[]>(c, grad);
	}
}
