package kip.javasim.molecular;

import static java.lang.Math.max;
import static java.lang.Math.sqrt;
import static scikit.numerics.Math2.sqr;

import java.util.ArrayList;

import scikit.util.Utilities;

public class BinnedPts3d<A extends Pt3d> {
	private double _L; // system length
	private int _cols; // number of bins along length of volume
	private double _dx; // distance between neighboring bins
	private ArrayList<A>[] _cells;
	
	
	@SuppressWarnings("unchecked")
	public BinnedPts3d(double L, int cols) {
		_L = L;
		_cols = max(cols, 1);
		_dx = L / _cols;
		_cells = new ArrayList[_cols*_cols*_cols];
		for (int i = 0; i < _cols*_cols*_cols; i++)
			_cells[i] = new ArrayList<A>();
	}
	
	public void addAtom(A p) {
		if (!isValidPoint(p.x, p.y, p.z))
			throw new IllegalArgumentException("Atom out of bounds");
		_cells[pointToIndex(p.x, p.y, p.z)].add(p);
	}
	
	public void removeAtom(A p) {
		if (!isValidPoint(p.x, p.y, p.z))
			throw new IllegalArgumentException("Atom out of bounds");
		_cells[pointToIndex(p.x, p.y, p.z)].remove(p);
	}
	
	public ArrayList<A> getAll() {
		ArrayList<A> ret = new ArrayList<A>();
		for (ArrayList<A> cell : _cells)
			ret.addAll(cell);
		return ret;
	}
	
	ArrayList<A> tempArray = new ArrayList<A>();
	
	public ArrayList<A> atomsWithinRange(double x, double y, double z, double R) {
		// imax represents the maximum number of simulation cells that _R can span 
		int imax = (int)(R/_dx+1.0);
		// make sure that the atoms in cells are not double counted
		if (2*imax+1 > _cols)
			return atomsWithinRangeSlow(x, y, z, R);
		
		x = (x+_L)%_L;
		y = (y+_L)%_L;
		z = (z+_L)%_L;
		int i1 = coordToIndex(x);
		int j1 = coordToIndex(y);
		int k1 = coordToIndex(z);
		
		tempArray.clear();
		ArrayList<A> ret = tempArray;
//		ArrayList<A> ret = new ArrayList<A>();
		
		// if the center-to-center distance squared of two cells exceeds d2Cutoff then no atoms within the
		// cells can interact.
		// sqrt(3.) represents the corner-to-corner distance of a single cell
		int d2Cutoff = (int) (sqr(R/_dx+sqrt(3.))+1e-8);
		
		for (int di = -imax; di <= imax; di++) {
			for (int dj = -imax; dj <= imax; dj++) {
				for (int dk = -imax; dk <= imax; dk++) {
					if (di*di + dj*dj + dk*dk <= d2Cutoff) {
						int i2 = i1+di;
						int j2 = j1+dj;
						int k2 = k1+dk;
						i2 = (i2+_cols)%_cols;
						j2 = (j2+_cols)%_cols;						
						k2 = (k2+_cols)%_cols;

						// it is significantly faster to loop i explicitly, rather than
						// for (A p : cell)
						ArrayList<A> cell = _cells[_cols*_cols*k2+_cols*j2+i2];
						for (int i = 0; i < cell.size(); i++) {
							A p = cell.get(i);
							double dx = p.x - x + (i1+di-i2)*_dx;
							double dy = p.y - y + (j1+dj-j2)*_dx;
							double dz = p.z - z + (k1+dk-k2)*_dx;
							if (dx*dx + dy*dy + dz*dz < R*R)
								ret.add(p);
						}
					}
				}
			}
		}
		
		boolean TESTING = true;
		if (TESTING && ret.size() != atomsWithinRangeSlow(x, y, z, R).size())
			throw new IllegalStateException("Counting error.");
		return ret;
	}
	
	// for testing purposes only
	public ArrayList<A> atomsWithinRangeSlow(double x, double y, double z, double R) {
		ArrayList<A> ret = new ArrayList<A>();
		
		for (ArrayList<A> cell : _cells) {
			for (int i = 0; i < cell.size(); i++) {
				A p = cell.get(i);
				double dx = x - p.x;
				double dy = y - p.y;
				double dz = z - p.z;
				dx = Utilities.periodicOffset(_L, dx);
				dy = Utilities.periodicOffset(_L, dy);
				dz = Utilities.periodicOffset(_L, dz);
				if (dx*dx + dy*dy + dz*dz < R*R)
					ret.add(p);
			}
		}
		return ret;
	}
	
	
	private boolean isValidPoint(double x, double y, double z) {
		return ((0 <= x && x < _L) && (0 <= y && y < _L) && (0 <= z && z < _L));
	}

	// rounding errors here are OK, as long as they occur in just this
	// one function	
	private int coordToIndex(double x) {
		int i = (int)(x/_dx);
		assert (0 <= i && i < _cols);
		return i;
	}
	
	private int pointToIndex(double x, double y, double z) {
		return coordToIndex(z)*_cols*_cols + coordToIndex(y)*_cols + coordToIndex(x);
	}
	
	// returns center of grid element
//	private Tuple3d indexToPoint(int index) {
//		int i = index%_cols;
//		int j = (index/_cols)%_cols;
//		int k = index/(_cols*_cols);
//		return new Tuple3d((i+0.5)*_dx, (j+0.5)*_dx, (k+0.5)*_dx);
//	}	


}
