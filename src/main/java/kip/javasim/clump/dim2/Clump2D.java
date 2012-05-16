package kip.javasim.clump.dim2;

import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.exp;
import static java.lang.Math.floor;
import static java.lang.Math.min;
import static scikit.numerics.Math2.hypot;
import scikit.dataset.Accumulator;
import scikit.jobs.Job;
import scikit.jobs.params.Parameters;
import scikit.numerics.fft.FFT2D;


public class Clump2D extends AbstractClump2D {	
	PtsGrid pts;
	int t_cnt, numPts;
	double dt;
	double[] ptsX, ptsY;
	FFT2D fft;
	
	public Clump2D(Parameters params) {
		random.setSeed(params.iget("Random seed", 0));

		R = params.fget("R");
		L = params.fget("L");
		T = params.fget("T");
		dx = params.fget("dx");
		dt = params.fget("dt");
		
		numPts = (int)(L*L);
		pts = new PtsGrid(L, R, dx);
		ptsX = new double[numPts];
		ptsY = new double[numPts];
		randomizePts();
		t_cnt = 0;
		
		fft = new FFT2D((int)(2*L), (int)(2*L));
		fft.setLengths(L, L);
	}
	
	private void randomizePts() {
		for (int i = 0; i < numPts; i++) {
			ptsX[i] = random.nextDouble()*L;
			ptsY[i] = random.nextDouble()*L;			
			pts.add(ptsX[i], ptsY[i]);
		}
	}
	
	public void readParams(Parameters params) {
		T = params.fget("T");
		dt = params.fget("dt");
	}
	
	double dist2(double dx, double dy) {
		dx = abs(dx);
		dx = min(dx, L-dx);
		dy = abs(dy);
		dy = min(dy, L-dy);
		return dx*dx + dy*dy;
	}
	
	int slowCount(double x, double y) {
		int acc = 0;
		for (int i = 0; i < numPts; i++) {
			if (dist2(ptsX[i]-x, ptsY[i]-y) < R*R) {
				acc++;
			}
		}
		return acc;
	}
	
	void mcsTrial() {
		int i = random.nextInt(numPts);
		double x = ptsX[i];
		double y = ptsY[i];
		
		double xp = x + R*(2*random.nextDouble()-1);
		double yp = y + R*(2*random.nextDouble()-1);
		xp = (xp+L)%L;
		yp = (yp+L)%L;
//		assert(pts.countOverlaps(xp,yp) == slowCount(xp,yp));
//		assert(pts.countOverlaps(x,y) == slowCount(x,y));
		
		double dE = (pts.countOverlaps(xp,yp)-pts.countOverlaps(x,y))/(PI*R*R);
		if (dE < 0 || random.nextDouble() < exp(-dE/T)) {
			ptsX[i] = xp;
			ptsY[i] = yp;			
			pts.remove(x, y);
			pts.add(xp, yp);
		}
		t_cnt++;
	}
	
	public void simulate() {
		for (int i = 0; i < numPts*dt; i++) {
			mcsTrial();
			Job.yield();
		}
	}
	
	public Accumulator newStructureAccumulator(double binWidth) {
		// round binwidth down so that it divides KR_SP without remainder.
		binWidth = KR_SP / floor(KR_SP/binWidth);
		return new Accumulator(binWidth);
	}
	
	public void accumulateStructure(final Accumulator sf) {
		double[] scratch = fft.getScratch();
		int Lp = fft.dim1;
		performCoarseGraining(scratch, Lp);
		fft.transform(scratch, new FFT2D.MapFn() {
			public void apply(double k1, double k2, double re, double im) {
				double k = hypot(k1, k2);
				double kR = k*R;
				if (kR > 0 && kR <= 4*KR_SP)
					sf.accum(kR, (re*re+im*im)/(L*L));
			}
		});
	}	
	
	public double[] coarseGrained() {
		return pts.rawElements;
	}
	
	public int numColumns() {
		return pts.gridCols;
	}
	
	public double time() {
		return (double)t_cnt/numPts;
	}
	
	private void performCoarseGraining(double[] field, int Lp) {
		double dx = L/Lp;
		for (int i = 0; i < Lp*Lp; i++)
			field[i] = 0;
		for (int k = 0; k < ptsX.length; k++) {
			int i = (int)(Lp*ptsX[k]/L);
			int j = (int)(Lp*ptsY[k]/L);
			assert(i < Lp && j < Lp);
			field[Lp*j+i] += 1/(dx*dx);
		}
	}
}
