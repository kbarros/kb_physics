package kip.javasim.clump.dim3;
import static java.lang.Math.*; 


class PtsList3D {
	double L, R;
	double[] xs;
	double[] ys;
	double[] zs;
	int cnt;
	static final int CELL_SIZE = 8;
	
	public PtsList3D(double L, double R) {
		this.L = L;
		this.R = R;
		
		xs = new double[CELL_SIZE];
		ys = new double[CELL_SIZE];		
		zs = new double[CELL_SIZE];
		cnt = 0;
	}
	
	private void resize(int n) {
		double[] nxs = new double[n];
		double[] nys = new double[n];
		double[] nzs = new double[n];
		System.arraycopy(xs, 0, nxs, 0, cnt);		
		System.arraycopy(ys, 0, nys, 0, cnt);
		System.arraycopy(zs, 0, nzs, 0, cnt);
		xs = nxs;
		ys = nys;
		zs = nzs;
	}
	

	
	public void add(double x, double y, double z) {
		if (cnt == xs.length)
			resize(2*xs.length);
		xs[cnt] = x;
		ys[cnt] = y;
		zs[cnt] = z;
		cnt++;
	}
	
	public void remove(double x, double y, double z) {
		if (xs.length > CELL_SIZE && cnt < xs.length/4)
			resize(xs.length/2);
		for (int i = 0; i < cnt; i++) {
			if (xs[i] == x && ys[i] == y && zs[i] == z) {
				xs[i] = xs[cnt-1];
				ys[i] = ys[cnt-1];
				zs[i] = zs[cnt-1];
				cnt--;
				return;
			}
		}
		throw new IllegalArgumentException();
	}

	public int countOverlaps(double x, double y, double z) {
		if (cnt == 0) return 0;
		double dx = x-xs[0];
		double dy = y-ys[0];
		double dz = z-zs[0];
		if (dx >  L/2) x -= L;
		if (dx < -L/2) x += L;
		if (dy >  L/2) y -= L;
		if (dy < -L/2) y += L;
		if (dz >  L/2) z -= L;
		if (dz < -L/2) z += L;
		int acc = 0;
		for (int i = 0; i < cnt; i++) {
			dx = x-xs[i];
			dy = y-ys[i];
			dz = z-zs[i];
			if (dx*dx + dy*dy + dz*dz < R*R)
				acc++;
		}
		return acc;
	}
}


public class PtsGrid3D {
	double L; // system length
	double R; // interaction range 
	
	int gridCols;  // columns in grid  
	double dx; // distance between grid elements
	PtsList3D[] grid; // grid of cells, each containing list of particle positions
	double[] rawElements; // particle density stored for each cell
	
	public PtsGrid3D(double L, double R, double dx) {
		this.L = L;
		this.R = R;
		this.dx = dx;
		
		gridCols = (int)(L/dx);
		this.dx = L / gridCols;
		assert (gridCols > 3);
		assert (4*R <= L);
		
		// first grid element has  (x in [0, dx),   y,z in [0, dx))
		// second grid element has (x in [dx, 2*dx) y,z in [0, dx))
		// ...
		// ...                     (x,z in [0, dx), y in [dx, 2*dx])
		// ...
		grid = new PtsList3D[gridCols*gridCols*gridCols];
		for (int i = 0; i < grid.length; i++)
			grid[i] = new PtsList3D(L, R);
		rawElements = new double[gridCols*gridCols*gridCols];
	}
	
	private int gridIndex(double x, double y, double z) {
		int i = (int)(x/dx);
		int j = (int)(y/dx);
		int k = (int)(z/dx);
		assert(i < gridCols && j < gridCols && k < gridCols);
		return k*gridCols*gridCols + j*gridCols + i;
	}
	
	public int countOverlaps(double x, double y, double z) {
		// cubic cells are represented by point at the center (cx,cy,cz)
		// distance from center to corner is dx*sqrt(3)/2
		// if the center of a cell is within distance R+dx*sqrt(3)/2 from a point
		// then it might contain overlapping particles.
		// the extra epsilon is to deal with floating point rounding.
		double r_lo = R-dx*sqrt(3)/2-1e-8;
		double r_hi = R+dx*sqrt(3)/2+1e-8;
		
		// center of cell which contains (x,y,z)
		double c1x = dx*(floor(x/dx) + 0.5);
		double c1y = dx*(floor(y/dx) + 0.5);
		double c1z = dx*(floor(z/dx) + 0.5);
		
		int acc = 0;
		
		int imax = (int)ceil(R/dx);
		for (int i = -imax; i <= imax; i++) {
			for (int j = -imax; j <= imax; j++) {
				for (int k = -imax; k <= imax; k++) {
					// center of cell to test
					double c2x = c1x + dx*i;
					double c2y = c1y + dx*j;
					double c2z = c1z + dx*k;
					double d2 = (x-c2x)*(x-c2x) + (y-c2y)*(y-c2y) + (z-c2z)*(z-c2z);
					
					if (d2 < r_hi*r_hi) {
						int idx = gridIndex((c2x+L)%L, (c2y+L)%L, (c2z+L)%L);
						if (d2 < r_lo*r_lo)
							acc += grid[idx].cnt;
						else
							acc += grid[idx].countOverlaps(x,y,z);
					}
				}
			}
		}
		return acc;
	}
	
	public void remove(double x, double y, double z) {
		grid[gridIndex(x, y, z)].remove(x, y, z);
		rawElements[gridIndex(x,y,z)] -= 1/(dx*dx*dx);
	}
	
	public void add(double x, double y, double z) {
		grid[gridIndex(x, y, z)].add(x, y, z);
		rawElements[gridIndex(x,y,z)] += 1/(dx*dx*dx);
	}
}