package kip.javasim.fun;

import static java.lang.Math.*;

import java.awt.Color;

import kip.javasim.Random;
import scikit.dataset.Function;
import scikit.dataset.Histogram;
import scikit.dataset.PointSet;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;



class LinearRegression {
	public double a;
	public double b;
	
	// Fit data with expression of the form "y = ax + b"
	public LinearRegression(double[] x, double[] y, double[] error, double xlo, double xhi) {
		int n = 0;
		double mx=0, my=0, mxx=0, mxy=0;
		
		for (int i = 0; i < x.length; i++) {
			if (xlo <= x[i] && x[i] < xhi) {
				mx += x[i];
				my += y[i];
				mxx += x[i]*x[i];
				mxy += x[i]*y[i];
				n += 1;
			}
		}
		
		if (n < 2) {
			a = b = 0;
		}
		else {
			mx /= n;
			my /= n;
			mxx /= n;
			mxy /= n;
			a = (mxy - mx*my) / (mxx - mx*mx);
			b = (mxx*my - mxy*mx) / (mxx - mx*mx);
		}
	}
}


public class WindingWalkApp extends Simulation {
	Plot escapePlot = new Plot("Fields");
	Plot rangePlot = new Plot("Equilibrium distribution of x=theta/alpha");
	Histogram escapeAccum;
	Histogram rangeAccum;
	Random random;
	
	public static void main(String[] args) {
		new Control(new WindingWalkApp(), "WindingWalk");
	}


	public void load(Control c) {
		c.frame(escapePlot);
		c.frame(rangePlot);
		params.add("Max time", "1e8");
		params.add("Bin width", 0.01);
		params.add("Dynamics", new ChoiceValue(new String[]{"Fast", "Slow"}));
        params.add("Alpha", 1.5*Math.PI);
		params.add("dt", 0.01);
		params.add("Random seed", 0);
		params.addm("Min fit cutoff", "1e5");
		params.addm("Max fit cutoff", "1e8");
		
		params.add("Survival exponent");
	}
	
	public void animate() {
		//
		// extract slope from log-log distribution of first passage times 
		//
		double[][] data = escapeAccum.copyData().columns();
		double[] x = data[0];
		double[] y = data[1];
		for (int i = 0; i < y.length; i++) {
			// x[i] = pow(10, x[i]);
			double scaledRealTimeBinWidth = pow(10, x[i]);  
			y[i] = log10(y[i]/scaledRealTimeBinWidth);
		}
		PointSet logEscapeAccum = new PointSet(x, y);
		
		double loCutoff = params.fget("Min fit cutoff");
		double hiCutoff = params.fget("Max fit cutoff");
		final LinearRegression lr = new LinearRegression(x, y, null, log10(loCutoff), log10(hiCutoff));
		Function f = new Function() {
			public double eval(double x) {
				return lr.a*x + lr.b;
			}
		};
		params.set("Survival exponent", lr.a);
		
		escapePlot.registerLines("Escape times", logEscapeAccum, Color.RED);
		escapePlot.registerLines("Regression", f, Color.BLUE);
		
		//
		// plot distribution of surviving angles
		//
		rangePlot.registerLines("Angle distribution", rangeAccum, Color.RED);
		rangePlot.registerLines("Theoretical", new Function(-1, 1) {
			public double eval(double x) { return (PI/4)*cos(PI*x/2); }
		}, Color.BLUE);
	}
    
	public void clear() {
		escapePlot.clear();
		rangePlot.clear();
	}
	
	double randomDir(Random rand) {
		return 2.0*rand.nextInt(2) - 1.0;
	}
	
	double sqr(double x) {
		return x*x;
	}

	public void run() {
		String dynamics = params.sget("Dynamics");
		boolean slow = dynamics.equals("Slow");
		double alpha = params.fget("Alpha");
		double initDt = params.fget("dt");
		double maxTime = params.fget("Max time");
		
		double binWidth = params.fget("Bin width");
		long randomSeed = params.iget("Random seed");
		Random rand = new Random(randomSeed);
		
		double diskRadius = 1.0;
		escapeAccum = new Histogram(binWidth);
		rangeAccum = new Histogram(binWidth);
		rangeAccum.setNormalizing(true);
		
		while (true) {
			// begin random walk
			
//			double r = 1*diskRadius;
//			double theta = 0;
//			double time = 0;
//			while (Math.abs(theta) < alpha && time < maxTime) {
//				double dt = initDt * (slow ? 1 : Math.max(1, sqr(r/diskRadius-1)));
//				theta += (Math.sqrt(dt) / r) * randomDir(rand);
//				r += (Math.sqrt(dt)) * randomDir(rand);
//				r = Math.max(r, diskRadius);
//				time += dt;
//			}

			double x = 10*diskRadius;
			double y = 0;
			double time = 0;
			int wrap = 0;
			double theta = atan2(y, x) + 2*PI*wrap;

			while (abs(theta) < alpha && time < maxTime) {
				double dt = initDt * (slow ? 1 : (x*x + y*y));
				double mag = sqrt(dt);
				double yp = y;
				
				switch (rand.nextInt(4)) {
				case 0: x += mag; break;
				case 1: x -= mag; break;
				case 2: y += mag; break;
				case 3: y -= mag; break;
				}
				
				// this is a singularity we're not supposed to hit, but it's possible with lattice dynamics;
				// interpret this point as infinite winding
				if (x == 0 && y == 0) {
					//assert(false);
					theta = Double.POSITIVE_INFINITY;
					break;
				}
				
				if (x <= 0) {
					// atan2 branches at (x < 0, y = +-0); make sure that
					// (y >= 0) is non-negative and (y < 0) is negative
					if (y == -0)
						y = +0;
					
					if (yp >= 0 && y < 0) // cross from positive to negative
						wrap++;
					if (yp < 0 && y >= 0) // cross from negative to positive
						wrap--;
				}
				time += dt;
				theta = atan2(y, x) + 2*PI*wrap;
			}

			
			if (binWidth < time && time < maxTime)
				escapeAccum.accum(log10(time));
			if (time >= maxTime)
				rangeAccum.accum(theta/alpha);
			
			Job.animate();
		}
	}
}