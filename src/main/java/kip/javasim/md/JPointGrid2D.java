package kip.javasim.md;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.sqrt;
import static scikit.numerics.Math2.sqr;

import java.util.ArrayList;

import scikit.util.Point;
import scikit.util.Utilities;

public class JPointGrid2D<Pt extends Point> {
	private double _L;
	private int _cols;
	private double _dx;
	private boolean _periodic;
	private Pt[] _points;
	private ArrayList<Pt>[] _cells;
	
	
	@SuppressWarnings("unchecked")
	public JPointGrid2D(double L, int cols, boolean periodic, Pt[] points) {
		_L = L;
		_cols = max(cols, 1);
		_dx = L / _cols;
		_periodic = periodic;
		_points = points;
		_cells = new ArrayList[_cols*_cols];
		for (int i = 0; i < _cols*_cols; i++)
			_cells[i] = new ArrayList<Pt>();
		initialize();
	}

	public void initialize() {
		for (ArrayList<Pt> c : _cells)
			c.clear();
		for (Pt p : _points) {
			_cells[pointToIndex(p.x, p.y)].add(p);
		}
	}
	
	ArrayList<Pt> tempArray = new ArrayList<Pt>();
	
	public ArrayList<Pt> pointOffsetsWithinRange(Point p, double R) {
		int imax = (int)(R/_dx+1.0);
		if (2*imax+1 > _cols)
			return pointOffsetsWithinRangeSlow(p, R);
		
		double x = (p.x+_L)%_L;
		double y = (p.y+_L)%_L;
		int index = pointToIndex(x, y);
		int i1 = index%_cols;
		int j1 = index/_cols;
		
		tempArray.clear();
		ArrayList<Pt> ret = tempArray;
//		ArrayList<Pt> ret = new ArrayList<Pt>();
		
		int d2Cutoff = (int) (sqr(R/_dx+sqrt(2))+1e-8);
		
		for (int di = -imax; di <= imax; di++) {
			for (int dj = -imax; dj <= imax; dj++) {
				if (di*di + dj*dj <= d2Cutoff) {
					int i2 = i1+di;
					int j2 = j1+dj;
					if (_periodic) {
						i2 = (i2+_cols)%_cols;
						j2 = (j2+_cols)%_cols;						
					}
					else if (min(i2,j2) < 0 || max(i2,j2) >= _cols) {
						continue;
					}
					
					// it is significantly faster to loop i explicitly, rather than
					// for (Pt p : cell)
					ArrayList<Pt> cell = _cells[_cols*j2+i2];
					for (int i = 0; i < cell.size(); i++) {
						Pt p2 = cell.get(i);
						double dx = p2.x - p.x + (i1+di-i2)*_dx;
						double dy = p2.y - p.y + (j1+dj-j2)*_dx;
						if (dx*dx + dy*dy < R*R)
							ret.add(p2);
					}
				}
			}
		}
//		if (ret.size() != pointOffsetsWithinRangeSlow(p, R).size())
//			throw new IllegalStateException("Counting error.");
		return ret;
	}
	
	// for testing purposes only
	public ArrayList<Pt> pointOffsetsWithinRangeSlow(Point p, double R) {
		ArrayList<Pt> ret = new ArrayList<Pt>();
		
		for (int i = 0; i < _points.length; i++) {
			double dx = p.x - _points[i].x;
			double dy = p.y - _points[i].y;
			if (_periodic) {
				dx = Utilities.periodicOffset(_L, dx);
				dy = Utilities.periodicOffset(_L, dy);
			}
			if (dx*dx + dy*dy < R*R)
				ret.add(_points[i]);
		}
		return ret;
	}
	
	
	// rounding errors here are OK, as long as they occur in just this
	// one function
	private int pointToIndex(double x, double y) {
		int i = (int)(x/_dx);
		int j = (int)(y/_dx);
		assert(i < _cols && j < _cols);
		return j*_cols+i;
	}
	
	// returns center of grid element
	@SuppressWarnings("unused")
	private Point indexToPoint(int index) {
		int i = index%_cols;
		int j = index/_cols;
		return new Point((i+0.5)*_dx, (j+0.5)*_dx);
	}	
}
