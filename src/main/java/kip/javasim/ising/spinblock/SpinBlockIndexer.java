package kip.javasim.ising.spinblock;

import static java.lang.Integer.*;


/**
class SpinBlockIndexer

The heart of this class is the fillArray() method.  fillArray() takes a position,
x, and returns an array of "spin-block indices" which together account for all spins
in the range, R, of x.  Each spin-block index is actually a pair of numbers: a scale,
and a block index at that scale.  The finest scale, 0, is at the level of spins
themselves.  The coarsest scale is given by maxScale, and is roughly the size of the
range itself.

This class uses state internally, and is therefore not thread safe.  Nonetheless, the
semantics of it's methods should not depend on previous state.
*/

class SpinBlockIndexer {
	private int L;
	private int maxScale;
	
	private int[] indices;
	private int numIndices;
	
	
	public SpinBlockIndexer(int L, int R) {
		assert (L == highestOneBit(L));
		assert (L >= 2*R+1);
		this.L = L;
		maxScale = numberOfTrailingZeros(2*highestOneBit(2*R+1));
	}
	
	
	public int maxScale() {
		return maxScale;
	}
	
	
	public int[] newArray() {
		return new int[4*maxScale+1];
	}
	
	
	private void fillArrayAux2(int xlo, int xhi, int scale) {
		assert (scale >= 0);
		assert (xlo >> scale == xhi >> scale);
		
		if (xhi+1 - xlo == 1 << scale) {
			indices[numIndices++] = scale;
			indices[numIndices++] = xlo >> scale;
			return;
		}
		
		int scale_m = scale - 1;
		int xlo_scaled = xlo >> scale_m;
		int xhi_scaled = xhi >> scale_m;
		
		if (xlo_scaled == xhi_scaled) {
			fillArrayAux2(xlo, xhi, scale_m);
		}
		else {
			assert (xlo_scaled + 1 == xhi_scaled);
			int xmid = xhi_scaled << scale_m;
			fillArrayAux2(xlo, xmid-1, scale_m);
			fillArrayAux2(xmid, xhi, scale_m);
		}
	}
	
	
	private void fillArrayAux1(int xlo, int xhi) {
		int xlo_scaled = xlo >> maxScale;
		int xhi_scaled = xhi >> maxScale;
		
		if (xlo < 0) {
			fillArrayAux1(0, xhi);
			fillArrayAux1((xlo+L)%L, L-1);
		}
		else if (xhi >= L) {
			fillArrayAux1(xlo, L-1);
			fillArrayAux1(0, (xhi+L)%L);
		}
		else if (xlo_scaled < xhi_scaled) {
			assert (xlo_scaled + 1 == xhi_scaled);
			int xmid = xhi_scaled << maxScale;
			fillArrayAux1(xlo, xmid-1);
			fillArrayAux1(xmid, xhi);
		}
		else {
			fillArrayAux2(xlo, xhi, maxScale);
		}
	}
	
	
	public void fillArray(int xlo, int xhi, int[] result) {
		indices = result;
		numIndices = 0;
		fillArrayAux1(xlo, xhi);
		indices[numIndices++] = -1;
	}
}
