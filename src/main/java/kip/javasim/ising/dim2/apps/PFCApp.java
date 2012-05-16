package kip.javasim.ising.dim2.apps;

import static scikit.numerics.Math2.cube;
import static scikit.numerics.Math2.sqr;
import static scikit.util.Utilities.format;
import kip.javasim.Random;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DoubleValue;
import scikit.jobs.params.Parameters;
import scikit.numerics.fft.FFT2D;
import scikit.numerics.fn.Function2D;

class PFC2D {
	Random random = new Random();
	double L, dx;
	int Lp;
	double dt, t; 
	double r, da, phi0;
	double alpha, beta, gamma;
	double phibar, amplitude;
	double tsnap, strainRate;

	boolean typA = false;
	boolean typB = true;
	int ntyps = (typA ? 1 : 0) + (typB ? 1 : 0);
	
	double boundaryFrac = 0.05;
	double strainMag = 0.5;

	double[] rho, rhoPrev, conc, concPrev;
	double[] rhoSnap, concSnap;
	double[] scratch1, scratch2, scratch3;
	
	FFT2D fft;
	
	public PFC2D(Parameters params) {
		random.setSeed(params.iget("Random seed", 0));
		
		L = params.fget("L");
		dx = params.fget("dx");
		Lp = (int) (L / dx);
		dx = L / Lp;
		params.set("dx", dx);
		
		readParams(params);
		
		rho = new double[Lp*Lp];
		rhoPrev = new double[Lp*Lp];
		conc = new double[Lp*Lp];
		rhoSnap = new double[Lp*Lp];
		concSnap = new double[Lp*Lp];
		scratch1 = new double[Lp*Lp];
		scratch2 = new double[Lp*Lp];
		scratch3 = new double[Lp*Lp];
		fft = new FFT2D(Lp, Lp);
		fft.setLengths(L, L);
		
		t = 0;
		tsnap = 0;
		
		phibar = params.fget("Initial density");
		amplitude = (4./5.)*(phibar + Math.sqrt(-15*r - 36*phibar*phibar));
		System.out.println(amplitude + phi0);
		
		// buildRandomInitial();
		buildLayeredInitial();
		buildSnapshot();
	}
	
	public void readParams(Parameters params) {
		r = params.fget("r");
		dt = params.fget("dt");
		
		double g = params.fget("a1/a0");
		da = (g*g-1)/(g*g+1); 
		
		phi0 = params.fget("phi0");
		alpha = params.fget("alpha");
		beta = params.fget("beta");
		gamma = params.fget("gamma");
		strainRate = params.fget("Strain rate");
	}
	
	public void buildRandomInitial() {
		for (int i = 0; i < Lp*Lp; i++) {
			rho[i] = phibar + phi0 + 0.1*random.nextGaussian();
			rhoPrev[i] = rho[i];
			
			int xp = i % Lp;			
			if (xp > (1-2*boundaryFrac)*Lp)
				conc[i] = 1;
			else if (xp < (0+2*boundaryFrac)*Lp)
				conc[i] = -1;
			else
				conc[i] = 0.1*random.nextGaussian();
		}
	}
	
	public void buildLayeredInitial() {
		for (int i = 0; i < Lp*Lp; i++) {
			rho[i] = phibar + phi0 + 0.1*random.nextGaussian();
			rhoPrev[i] = rho[i];
			
			int xp = i % Lp;
			int yp = i / Lp;
			if (xp > (1-2*boundaryFrac)*Lp)
				conc[i] = 1;
			else if (xp < (0+2*boundaryFrac)*Lp)
				conc[i] = -1;
			else {
				conc[i] = 0.1*random.nextGaussian() + (yp > Lp/2 ? 1 : -1) ;
			}
		}
	}

	public void buildSnapshot() {
		for (int i = 0; i < Lp*Lp; i++) {
			int xp = i % Lp;
			int yp = i / Lp;
			
			double fudgeScale = 1.2;
			if (xp > (1-boundaryFrac)*Lp) {
				double sign = +1;
				double R = fudgeScale*Math.sqrt(1 + sign*da);
				int n = (int)Math.round(L / (2*Math.PI*R));
				R = L / (2*Math.PI*n);
				rhoSnap[i] = amplitude*Math.sin((yp*dx)/R) + phi0;
				
				concSnap[i] = sign;
			}
			else if (xp < (0+boundaryFrac)*Lp) {
				double sign = -1;
				double R = fudgeScale*Math.sqrt(1 + sign*da);
				int n = (int)Math.round(L / (2*Math.PI*R));
				R = L / (2*Math.PI*n);
				rhoSnap[i] = amplitude*Math.sin((yp*dx)/R) + phi0;
				concSnap[i] = sign;
			}
			else {
				concSnap[i] = rhoSnap[i] = 0;
			}
		}
		tsnap = t;
	}
	
	public void takeSnapshot() {
		for (int i = 0; i < Lp*Lp; i++) {
			rhoSnap[i] = rho[i];
			concSnap[i] = conc[i];
		}
		tsnap = t;
	}
	
	public void simulate() {
		stepDensity();
		stepConcentration();
		t += dt;
	}
	
	public void laplacian(double[] src, double[] dst) {
		fft.convolve(src, dst, new Function2D() {
			public double eval(double kx, double ky) {
				double k2 = kx*kx + ky*ky;
				return -k2;
			}
		});
	}
	
	public void regularizedLaplacian(double[] src, double[] dst) {
		fft.convolve(src, dst, new Function2D() {
			public double eval(double kx, double ky) {
				double k2 = kx*kx + ky*ky;
				return -k2 / (1.0 + k2*k2*k2);
			}
		});
	}
	
	
	
	public double R2(int i) {
		return 1 + da*conc[i];
	}
	
	// [1 + R2 \del^2] src
	// cannot use same array for src and dst
	public void operatorA(double src[], double dst[]) {
		laplacian(src, dst);
		for (int i = 0; i < Lp*Lp; i++) {
			dst[i] = src[i] + R2(i) * dst[i];
		}
	}
	
	// [1 + \del^2 R2] src
	// cannot use same array for src and dst
	public void operatorB(double src[], double dst[]) {
		for (int i = 0; i < Lp*Lp; i++)
			dst[i] = R2(i) * src[i];
		laplacian(dst, dst);
		for (int i = 0; i < Lp*Lp; i++)
			dst[i] = src[i] + dst[i];
	}
	
	// bilinear interpolation
	// (x,y)=(0,0) is the center of the leftmost cell 
	public double interpolate(double[] data, double xp, double yp) {
		// ensure that x,y are in the range (0, Lp)
		xp = (xp%Lp + Lp) % Lp;
		yp = (yp%Lp + Lp) % Lp;
		int x0 = (int)xp;
		int x1 = (x0+1)%Lp;
		int y0 = (int)yp;
		int y1 = (y0+1)%Lp;
		
		double ax = xp % 1;
		double ay = yp % 1;
		return ((1-ax)*(1-ay)*data[y0*Lp+x0] + 
				  (ax)*(1-ay)*data[y0*Lp+x1] +
				(1-ax)*  (ay)*data[y1*Lp+x0] + 
				  (ax)*  (ay)*data[y1*Lp+x1]); 
	}
	
	public void strainForce(double[] src, double[] target, double[] dst) {
		double strain = (t - tsnap) * strainRate;
		for (int i = 0; i < Lp*Lp; i++) {
			int x = i % Lp;
			int y = i / Lp;
			if (x > (1-boundaryFrac)*Lp) {
				dst[i] += strainMag*(src[i] - interpolate(target, x, y + strain/dx));
			}
		}

		for (int i = 0; i < Lp*Lp; i++) {
			int x = i % Lp;
			int y = i / Lp;
			if (x < (0+boundaryFrac)*Lp) {
				dst[i] += strainMag*(src[i] - interpolate(target, x, y - strain/dx));
			}
		}
	}
	
	public void stepDensity() {
		// r rho + (rho-rho_0)^3 + beta rho^3 (c^2 - 1)^2
		for (int i = 0; i < Lp*Lp; i++) {
			scratch1[i] = r*rho[i] + cube(rho[i]-phi0); // + beta*cube(rho[i])*sqr(sqr(conc[i])-1);
		}
		
		if (typA) {
			// [1 + \del^2 R^2] [1 + R^2 \del^2] rho
			operatorA(rho, scratch2);
			operatorB(scratch2, scratch3);
			for (int i = 0; i < Lp*Lp; i++)
				scratch1[i] += scratch3[i] / ntyps;  
		}
		
		if (typB) {
			// [1 + R^2 \del^2] [1 + \del^2 R^2] rho
			operatorB(rho, scratch2);
			operatorA(scratch2, scratch3);
			for (int i = 0; i < Lp*Lp; i++)
				scratch1[i] += scratch3[i] / ntyps;  
		}

		strainForce(rho, rhoSnap, scratch1);
		
		// \del^2 / (1 + \del^6) [ ... ]
		regularizedLaplacian(scratch1, scratch1);

//		// update according to d rho/ dt = f
//		for (int i = 0; i < Lp*Lp; i++)
//			rho[i] += dt*scratch1[i];
		
//		// update according to: (gamma d^2 rho / dt^2 + d rho/ dt = f)
//		// where scratch1 stores f.  first derivative discretized forward in time.
//		for (int i = 0; i < Lp*Lp; i++) {
//			// rho_{t+dt}
//			double rhoNext = (1/(gamma + dt)) * (rho[i]*(2*gamma + dt) - gamma*rhoPrev[i] + dt*dt*scratch1[i]);
//			rhoPrev[i] = rho[i];
//			rho[i] = rhoNext;
//		}
		
		// update according to: (gamma d^2 rho / dt^2 + d rho/ dt = f)
		// where scratch1 stores f.  first derivative discretized central in time.
		for (int i = 0; i < Lp*Lp; i++) {
			// rho_{t+dt}
			double rhoNext = (1/(2*gamma + dt)) * (rho[i]*(4*gamma) + rhoPrev[i]*(dt-2*gamma) + 2*dt*dt*scratch1[i]);
			rhoPrev[i] = rho[i];
			rho[i] = rhoNext;
		}

	}
	
	public void stepConcentration() {
		// beta rho^4 c (c^2 - 1)
		for (int i = 0; i < Lp*Lp; i++) {
			scratch1[i] = beta*sqr(sqr(rho[i]))*conc[i]*(sqr(conc[i])-1);
		}
		
		if (typA) {
			// da (\del^2 rho) [1 + R^2 \del^2] rho
			laplacian(rho, scratch2);
			operatorA(rho, scratch3);
			for (int i = 0; i < Lp*Lp; i++)
				scratch1[i] += da * scratch2[i] * scratch3[i] / ntyps;
		}
		if (typB) {
			// da rho \del^2 [1 + \del^2 R^2] rho
			operatorB(rho, scratch2);
			laplacian(scratch2, scratch3);
			for (int i = 0; i < Lp*Lp; i++)
				scratch1[i] += da * rho[i] * scratch3[i] / ntyps;
		}
		
		// - alpha \del^2 c
		laplacian(conc, scratch2);
		for (int i = 0; i < Lp*Lp; i++)
			scratch1[i] += - alpha * scratch2[i];
		
		strainForce(conc, concSnap, scratch1);
		
		// \del^2 / (1 + \del^6) [ ... ]
		regularizedLaplacian(scratch1, scratch1);
		
		
		for (int i = 0; i < Lp*Lp; i++)
			conc[i] += dt*scratch1[i];
	}
	
	// E = \int dx^2 [ phi (r + (1+\del^2)^2) phi / 2 + phi^4 / 4 ]  
	public double energy() {
		// scratch1 stores (r + (1+\del^2)^2) phi
		fft.convolve(rho, scratch1, new Function2D() {
			public double eval(double kx, double ky) {
				double k2 = kx*kx + ky*ky;
				return r + (1-k2)*(1-k2);
			}
		});
		double ret = 0;
		for (int i = 0; i < Lp*Lp; i++) {			
			ret += rho[i]*scratch1[i]/2;      // phi (r + (1+\del^2)^2) phi / 2 
			ret += rho[i]*rho[i]*rho[i]*rho[i] / 4;	// phi^4 / 4
		}
		
		// shift energy so that it is zero in the ground state
		// ret += Lp*Lp/4; 
		
		return ret*dx*dx / (L*L);
	}
}



public class PFCApp extends Simulation {
	Grid grid1 = new Grid("Density");
	Grid grid2 = new Grid("Concentration");
	PFC2D sim;
	
	public static void main(String[] args) {
		new Control(new PFCApp(), "Phase field crystal model");
	}

	public void load(Control c) {
		c.frame(grid1);
		c.frame(grid2);
		params.add("L", 64.0);
		params.add("dx", 1.0);
		params.add("Initial density", -0.285);
		params.addm("r", new DoubleValue(-0.3, -1, 0).withSlider());
		params.addm("a1/a0", new DoubleValue(1.5, 0.5, 2).withSlider());
		params.addm("phi0", 0.7);
		params.addm("alpha", 1.0); // 0.05
		params.addm("beta", 2.0); // 5.0
		params.addm("gamma", 1.0); // alpha/beta in the notation of stefanovich
		params.addm("Strain rate", 0.0005);
		params.addm("dt", 1.0);
		params.add("Time");
		params.add("Energy");
		
		flags.add("Snap");
	}
	
	public void animate() {
		sim.readParams(params);
		
		grid1.setAutoScale();
		grid1.setDrawRange(true);
		grid1.setScale(0, 1.8);
		grid2.setScale(-1, 1);
	
		int Lp = sim.Lp;
		grid1.registerData(Lp, Lp, sim.rho);
		grid2.registerData(Lp, Lp, sim.conc);
		
		params.set("Time", format(sim.t));
		params.set("Energy", format(sim.energy()));
		
		if (flags.contains("Snap"))
			sim.takeSnapshot();
		flags.clear();
	}
	
	public void clear() {
		grid1.clear();
		grid2.clear();
	}
	
	public double getTime() {
		return sim.t;
	}
	
	public void run() {
		sim = new PFC2D(params);
		// sim.randomize();
		Job.animate();
		
		while (true) {
			for (int i = 0; i < 1; i++)
				sim.simulate();
			Job.animate();
		}
	}
}
