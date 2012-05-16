package kip.javasim.ising;

import java.util.HashSet;

public class PercolationSite2d extends NewmanZiff {
	int lx, ly;
	boolean openBoundary;
	
	public PercolationSite2d(int lx, int ly, boolean openBoundary) {
		super(lx*ly);
		this.lx = lx;
		this.ly = ly;
		this.openBoundary = openBoundary;
	}
	
	public void occupyAndBondSites(int a[], int id) {
		assert (a.length == lx*ly);
		for (int i = 0; i < lx*ly; i++)
			if (a[i] == id)
				occupyAndBondSite(i);
	}
	
	public void occupyAndBondSite(int i) {
		super.occupySite(i);
		int x = i % lx;
		int y = i / lx;			
		if (openBoundary) {
			if (y < ly-1)
				tryBond(i, i+lx, 0, 1);
			if (y > 0)
				tryBond(i, i-lx, 0, -1);
			if (x < lx-1)
				tryBond(i, i+1, 1, 0);
			if (x > 0)
				tryBond(i, i-1, -1, 0);
		}
		else {
			int yp = (y+1)%ly;
			int ym = (y-1+ly)%ly;
			int xp = (x+1)%lx;
			int xm = (x-1+lx)%lx;
			tryBond(i, yp*lx+x, 0, 1);
			tryBond(i, ym*lx+x, 0, -1);
			tryBond(i, y*lx+xp, 1, 0);
			tryBond(i, y*lx+xm, -1, 0);
		}
	}
	
	public void findHomologies() {
		if (openBoundary) {
			wrap_horizontal = wrap_vertical = wrap_diagonal = false;

			HashSet<Integer> left = new HashSet<Integer>();
			HashSet<Integer> right = new HashSet<Integer>();
			HashSet<Integer> bottom = new HashSet<Integer>();
			HashSet<Integer> top = new HashSet<Integer>();

			for (int y = 0; y < ly; y++) {
				left.add(clusterIndex(y*lx+0));
				right.add(clusterIndex(y*lx+(lx-1)));
			}
			for (int x = 0; x < lx; x++) {
				bottom.add(clusterIndex(0*lx+x));
				top.add(clusterIndex((ly-1)*lx+x));
			}

			left.retainAll(right);
			wrap_horizontal = !left.isEmpty();

			bottom.retainAll(top);
			wrap_vertical = !bottom.isEmpty();
		}
	}
	
	private void tryBond(int i, int j, int dx, int dy) {
		if (isOccupied(j))
			addBond(i, j, dx, dy);
	}
}
