package kip.javasim.clump.dim3;

import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.exp;
import static java.lang.Math.floor;
import static java.lang.Math.min;
import static scikit.numerics.Math2.hypot;
import scikit.dataset.Accumulator;
import scikit.jobs.Job;
import scikit.jobs.params.Parameters;
import scikit.numerics.fft.FFT3D;


public class Clump3D extends AbstractClump3D {
	PtsGrid3D pts;
	int t_cnt, numPts;
	double dt;
	double[] ptsX, ptsY, ptsZ;
	FFT3D fft;
	double[] fftScratch;
	
	public Clump3D(Parameters params) {
		random.setSeed(params.iget("Random seed", 0));

		R = params.fget("R");
		L = params.fget("L");
		T = params.fget("T");
		dx = params.fget("dx");
		dt = params.fget("dt");
		
		numPts = (int)(L*L*L);
		pts = new PtsGrid3D(L, R, dx);
		ptsX = new double[numPts];
		ptsY = new double[numPts];
		ptsZ = new double[numPts];
		randomizePts();
		t_cnt = 0;
		
		int dim = (int)(2*L);
		fft = FFT3D.create(dim, dim, dim);
		fft.setLengths(L, L, L);
		fftScratch = new double[dim*dim*dim];
	}
	
	private void randomizePts() {
		for (int i = 0; i < numPts; i++) {
			ptsX[i] = random.nextDouble()*L;
			ptsY[i] = random.nextDouble()*L;
			ptsZ[i] = random.nextDouble()*L;
			pts.add(ptsX[i], ptsY[i], ptsZ[i]);
		}
	}
	
	public void readParams(Parameters params) {
		T = params.fget("T");
		dt = params.fget("dt");
	}
	
	double dist2(double dx, double dy, double dz) {
		dx = abs(dx);
		dx = min(dx, L-dx);
		dy = abs(dy);
		dy = min(dy, L-dy);
		dz = abs(dz);
		dz = min(dz, L-dz);
		return dx*dx + dy*dy + dz*dz;
	}
	
	int slowCount(double x, double y, double z) {
		int acc = 0;
		for (int i = 0; i < numPts; i++) {
			if (dist2(ptsX[i]-x, ptsY[i]-y, ptsZ[i]-z) < R*R) {
				acc++;
			}
		}
		return acc;
	}
	
	void mcsTrial() {
		int i = random.nextInt(numPts);
		double x = ptsX[i];
		double y = ptsY[i];
		double z = ptsZ[i];
		
		double xp = x + R*(2*random.nextDouble()-1);
		double yp = y + R*(2*random.nextDouble()-1);
		double zp = z + R*(2*random.nextDouble()-1);
		
		xp = (xp+L)%L;
		yp = (yp+L)%L;
		zp = (zp+L)%L;
		
//		assert(pts.countOverlaps(xp,yp,zp) == slowCount(xp,yp,zp));
//		assert(pts.countOverlaps(x,y,z) == slowCount(x,y,z));
		
		double dE = (pts.countOverlaps(xp,yp,zp)-pts.countOverlaps(x,y,z))/(4*PI*R*R*R/3);
		if (dE < 0 || random.nextDouble() < exp(-dE/T)) {
			ptsX[i] = xp;
			ptsY[i] = yp;
			ptsZ[i] = zp;
			pts.remove(x, y, z);
			pts.add(xp, yp, zp);
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
		int Lp = fft.dim1;
		performCoarseGraining(fftScratch, Lp);
		fft.transform(fftScratch, new FFT3D.MapFn() {
			public void apply(double k1, double k2, double k3, double re, double im) {
				double k = hypot(k1, k2, k3);
				double kR = k*R;
				if (kR > 0 && kR <= 4*KR_SP)
					sf.accum(kR, (re*re+im*im)/(L*L*L));
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
		for (int i = 0; i < Lp*Lp*Lp; i++)
			field[i] = 0;
		for (int n = 0; n < ptsX.length; n++) {
			int i1 = (int)(Lp*ptsX[n]/L);
			int i2 = (int)(Lp*ptsY[n]/L);
			int i3 = (int)(Lp*ptsZ[n]/L);
			assert(i1 < Lp && i2 < Lp && i3 < Lp);
			field[Lp*Lp*i3+Lp*i2+i1] += 1/(dx*dx*dx);
		}
	}
}
