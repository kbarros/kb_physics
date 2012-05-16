package kip.javasim.clump.dim2;

import kip.javasim.LatticeNeighbors;
import static java.lang.Math.*; 


class PtsList {
	double L, R;
	double[] xs;
	double[] ys;
	int cnt;
	static final int CELL_SIZE = 8;
	
	public PtsList(double L, double R) {
		this.L = L;
		this.R = R;
		
		xs = new double[CELL_SIZE];
		ys = new double[CELL_SIZE];		
		cnt = 0;
	}
	
	private void resize(int n) {
		double[] nxs = new double[n];
		double[] nys = new double[n];
		System.arraycopy(xs, 0, nxs, 0, cnt);		
		System.arraycopy(ys, 0, nys, 0, cnt);
		xs = nxs;
		ys = nys;
	}
	

	
	public void add(double x, double y) {
		if (cnt == xs.length)
			resize(2*xs.length);
		xs[cnt] = x;
		ys[cnt] = y;
		cnt++;
	}
	
	public void remove(double x, double y) {
		if (xs.length > CELL_SIZE && cnt < xs.length/4)
			resize(xs.length/2);
		for (int i = 0; i < cnt; i++) {
			if (xs[i] == x && ys[i] == y) {
				xs[i] = xs[cnt-1];
				ys[i] = ys[cnt-1];
				cnt--;
				return;
			}
		}
		throw new IllegalArgumentException();
	}

//	private double dist2(double dx, double dy) {
//		if (dx < -L/2) dx = L + dx;
//		if (dx >  L/2) dx = L - dx;
//		if (dy < -L/2) dy = L + dy;
//		if (dy >  L/2) dy = L - dy;
//		return dx*dx + dy*dy;
//	}
//	public int countOverlaps(double x, double y) {
//		int acc = 0;
//		for (int i = 0; i < cnt; i++)
//			if (dist2(x-xs[i], y-ys[i]) < R*R)
//				acc++;
//		return acc;
//	}

	public int countOverlaps(double x, double y) {
		if (cnt == 0) return 0;
		double dx = x-xs[0];
		double dy = y-ys[0];
		if (dx >  L/2) x -= L;
		if (dx < -L/2) x += L;
		if (dy >  L/2) y -= L;
		if (dy < -L/2) y += L;
		int acc = 0;
		for (int i = 0; i < cnt; i++) {
			dx = x-xs[i];
			dy = y-ys[i];
			if (dx*dx + dy*dy < R*R)
				acc++;
		}
		return acc;
	}
}


public class PtsGrid {
	double L; // system length
	double R; // interaction range 
	
	int gridCols;  // columns in grid  
	double dx; // distance between grid elements
	PtsList[] grid;  // grid of cells, each containing list of particle positions
	double[] rawElements; // particle density stored for each cell
	
	LatticeNeighbors neigh1, neigh2;
	int[][] nlist1, nlist2;
	
	
	public PtsGrid(double L, double R, double dx) {
		this.L = L;
		this.R = R;
		this.dx = dx;
		
		gridCols = (int)(L/dx);
		this.dx = L / gridCols;
		
		// first grid element has  (x in [0, dx),   y in [0, dy))
		// second grid element has (x in [dx, 2*dx) y in [0, dy))
		// ...
		grid = new PtsList[gridCols*gridCols];
		for (int i = 0; i < grid.length; i++)
			grid[i] = new PtsList(L, R);
		rawElements = new double[gridCols*gridCols];
		
		// effectiveR is how many lattice spacings (length dx) away we need to
		// look to find relevant cells.  the sqrt(2) is necessary to to account for
		// the fact that the distance between cells should be counted as the distance
		// of the two closest corners.  the extra epsilon is to deal with floating point
		// rounding.
		double r_lo = R/dx-sqrt(2)-1e-8;
		double r_hi = R/dx+sqrt(2)+1e-8;
		neigh1 = new LatticeNeighbors(gridCols, gridCols, 0, r_lo, LatticeNeighbors.Type.PERIODIC);
		neigh2 = new LatticeNeighbors(gridCols, gridCols, r_lo, r_hi, LatticeNeighbors.Type.PERIODIC);
		nlist1 = new int[gridCols*gridCols][];
		nlist2 = new int[gridCols*gridCols][];
	}
	
	// rounding errors here are ok, as long as they occur in just this
	// one function
	private int gridIndex(double x, double y) {
		int i = (int)(x/dx);
		int j = (int)(y/dx);
		assert(i < gridCols && j < gridCols);
		return j*gridCols+i;
	}
	
	public int countOverlaps(double x, double y) {
		int j1 = gridIndex(x,y);
		if (nlist1[j1] == null) {
			nlist1[j1] = neigh1.get(j1);
			nlist2[j1] = neigh2.get(j1);
		}
		int acc = 0;
		for (int j2 : nlist1[j1])
			acc += grid[j2].cnt;
		for (int j2 : nlist2[j1])
			acc += grid[j2].countOverlaps(x, y);
		return acc;
	}
	
	public void remove(double x, double y) {
		grid[gridIndex(x, y)].remove(x, y);
		rawElements[gridIndex(x,y)] -= 1/(dx*dx);
	}
	
	public void add(double x, double y) {
		grid[gridIndex(x, y)].add(x, y);
		rawElements[gridIndex(x,y)] += 1/(dx*dx);
	}
}
