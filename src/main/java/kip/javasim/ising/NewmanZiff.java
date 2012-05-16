package kip.javasim.ising;

public class NewmanZiff {
	// number of sites
	public int n;
	
	public enum Topology {HORIZONTAL, VERTICAL, CROSS, POINT};
	
	protected boolean wrap_horizontal = false;
	protected boolean wrap_vertical = false;
	protected boolean wrap_diagonal = false;
	
	private final int EMPTY = Integer.MIN_VALUE;
	
	// if (parent[i] = EMPTY), then site is unoccupied, otherwise
	// if (parent[i] <  0), then -parent[i] is size of this cluster
	// if (parent[i] >= 0), then parent[i] is the parent of i
	private int parent[];
	
	// dx[i], dy[i]: displacement FROM site i TO parent[i] 
	private short dx[], dy[];
	
	public NewmanZiff(int n) {
		this.n = n;
		parent = new int[n]; 
		dx = new short[n];
		dy = new short[n];
		
		for (int i = 0; i < n; i++) {
			parent[i] = EMPTY; // root site
			dx[i] = dy[i] = 0;
		}
	}
	
	public void occupySite(int i) {
		if (parent[i] == EMPTY)
			parent[i] = -1;
	}
	
	public void addBond(int i, int j, int dx_i2j, int dy_i2j) {
		occupySite(i);
		occupySite(j);
		
		compressPath(i);
		compressPath(j);

		if (clusterIndex(i) == clusterIndex(j)) {
			// i and j are neighor sites.  if their displacements to their
			// (shared) root differs, then the cluster has wrapped
			boolean horiz = dx[i] != dx_i2j + dx[j];
			boolean vert  = dy[i] != dy_i2j + dy[j];
			if (horiz && vert)
				wrap_diagonal |= true;
			else {
				wrap_horizontal |= horiz;
				wrap_vertical |= vert;
			}
		}
		else {
			if (clusterSize(i) <= clusterSize(j))
				mergeRoots(i, j, dx_i2j, dy_i2j);
			else
				mergeRoots(j, i, -dx_i2j, -dy_i2j);
		}
	}
	
	public boolean isBonded(int i, int j) {
		return isOccupied(i) && isOccupied(j) && clusterIndex(i) == clusterIndex(j);
	}
	
	public boolean isOccupied(int i) {
		return parent[i] != EMPTY;
	}
	
	public int clusterSize(int i) {
		return isOccupied(i) ? -parent[clusterIndex(i)] : 0;
	}
	
	public int clusterIndex(int i) {
		return (parent[i] < 0) ? i : clusterIndex(parent[i]);
	}

	public boolean horizontalHomology() {
		return wrap_horizontal && !wrap_vertical && !wrap_diagonal;
	}
	
	public boolean verticalHomology() {
		return wrap_vertical && !wrap_horizontal && !wrap_diagonal;
	}
	
	public boolean pointHomology() {
		return !wrap_horizontal && !wrap_vertical && !wrap_diagonal;
	}
	
	public boolean crossHomology() {
		return
			(wrap_horizontal && wrap_vertical) ||
			(wrap_horizontal && wrap_diagonal) ||
			(wrap_vertical && wrap_diagonal);
	}
	
	public Topology getTopology() {
		if (horizontalHomology()) return Topology.HORIZONTAL;
		if (verticalHomology()) return Topology.VERTICAL;
		if (crossHomology()) return Topology.CROSS;
		if (pointHomology()) return Topology.POINT;
		return null;
	}
	
	private void mergeRoots(int i, int j, int dx_i2j, int dy_i2j) {
		int r_i = clusterIndex(i);
		int r_j = clusterIndex(j);
		assert (r_i != r_j);
		
		parent[r_j] += parent[r_i]; // r_j is root for i and j clusters 
		parent[r_i] = r_j;          // link r_i to r_j
		dx[r_i] = (short)(-dx[i] + dx_i2j + dx[j]); // distance from r_i to r_j
		dy[r_i] = (short)(-dy[i] + dy_i2j + dy[j]);
	}
	
	private void compressPath(int i) {
		if (parent[i] >= 0 && parent[parent[i]] >= 0) {
			compressPath(parent[i]);
			dx[i] += dx[parent[i]];
			dy[i] += dy[parent[i]];
			parent[i] = parent[parent[i]];
		}
		assert (i == clusterIndex(i) || parent[i] == clusterIndex(i));
	}
}
