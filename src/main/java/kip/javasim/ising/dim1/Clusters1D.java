package kip.javasim.ising.dim1;

import kip.javasim.ising.spinblock.SpinBlocks1D;


public class Clusters1D {
	// CHECK ME
	// static private final int NONE = Integer.MIN_VALUE;
	
	public SpinBlocks1D spins;
	public int L;                 // size of lattice
	public int R;				  // interaction range
	public double bondProbability;
	public int direction = 1;	  // spin direction for cluster formation
	
	public int[] numClusters;     // number of clusters of size s (n_s)
	public int[] randomSite;      // temporary storage for list of random sites
	
	// the parent[] array serves two purposes: It stores the cluster size
	// when the site is the (unique) root.  Otherwise, it stores the index of
	// the site's "parent." We can find the root from an arbitrary occupied site
	// by recursively following the parent array. This recursion terminates
	// when we encounter a negative value in the parent array, which indicates we have found
	// the unique cluster root.
	//    if (parent[s] >= 0)          parent[s] is the parent site index
	//    if (0 > parent[s] > NONE)    s is the root of size -parent[s]
	//    if (parent[s] == NONE)       site s is vacant
	private int[] parent;
	
	
	public Clusters1D(SpinBlocks1D spins, double p) {
		this.spins = spins;
		L = spins.L;
		R = spins.R;
		
		numClusters = new int[L+1];
		parent      = new int[L];
		randomSite  = new int[L];
		
		bondProbability = p;
	}
	
	
	public void setBondProbability(double p) {
		bondProbability = p;
	}
	
	
	public void buildClusters(boolean fast) {
		// initially all sites are empty, and there are no clusters
		for (int size = 0; size <= L; size++)
			numClusters[size] = 0;
		
		// make a singlet cluster for every spin aligned with "direction"
		for (int site = 0; site < L; site++) {
			if (spins.get(site) == direction) {
				// store the new cluster's size in parent[].  The negative sign
				// distinguishes newSite as a root, with a size value.  (Positive values
				// correspond to non-root sites with index pointers).
				parent[site] = -1;
				numClusters[1]++;
			}
		}
		
		for (int site = 0; site < L; site++) {
			if (spins.get(site) == direction) {				
				// Merge site with occupied neighbors, according to a bond probability.
				// 'root' is the index of the merged cluster root at each step
				int root = findRoot(site);
				
				// Find clusters the fast way...
				if (fast) {
					int sumOfNeighborsInRange = spins.sumInRange(site, site+R) - direction;
					int parallelSpins = (direction*sumOfNeighborsInRange + R) / 2;
					int n = binomialDeviate(parallelSpins, bondProbability);					
					
					selectRandomSitesInRange(site, n);
					for (int j = 0; j < n; j++) {
						root = mergeRoots(root, findRoot(randomSite[j]));
					}
				}
				// ... or the slow way
				else {
					for (int j = 1; j <= R; j++) {
						int neighbor = (site+j) % L;
						if (spins.get(neighbor) == direction && Math.random() < bondProbability) {
							root = mergeRoots(root, findRoot(neighbor));
						}
					}
				}
			}
		}
		
		sanityCheck();
	}


	// given a site index s, return the site index representing the root
	// of the cluster to which s belongs.
	private int findRoot(int s) {
		if (parent[s] < 0)
			return s; // i is the root site (with size -parent[s])
		else
			// first link parent[s] to the cluster's root to improve performance
			// (path compression.)  then return this value.
			return parent[s] = findRoot(parent[s]);
	}
	
	
	// merge two root sites into one.  this represents cluster merging.
	// use the heuristic that the root of the smaller cluster points to
	// the root of the larger cluster, in order to improve performance.
	// remember that parent[root] stores negative cluster size.
	private int mergeRoots(int r1, int r2) {
		// clusters are uniquely identified by their root sites.  if they
		// are the same, then the clusters are already merged, and we need
		// do nothing
		if (r1 == r2)
			return r1;
		// if r1 has smaller cluster size than r2, reverse (r1,r2) labels
		else if (-parent[r1] < -parent[r2])
			return mergeRoots(r2, r1);
		else /* (-parent[r1] > -parent[r2]) */ {
			// update the cluster count, and second cluster moment to account for the
			// loss of two small clusters and the gain of one big cluster
			numClusters[-parent[r1]]--;
			numClusters[-parent[r2]]--;
			numClusters[-parent[r1]-parent[r2]]++;
			// the cluster at r1 now includes the sites of old cluster at r2
			parent[r1] += parent[r2];
			// make r1 the new parent of r2
			parent[r2] = r1;
			// return the new root site r1
			return r1;
		}
	}
	
	
	private int binomialDeviate(int N, double p) {
		if (N < 25) {
			int n = 0;
			for (int i = 0; i < N; i++)
				if (Math.random() < p) n++;
			return n;
		}
		else { /* if(N*p < 1.0) { */
			double g = Math.exp(-N*p);
			int n = -1;
			double t = 1;
			do {
				n++;
				t *= Math.random();
			} while (t > g);

			return n;
		}
/*		else {
			assert false : "Can't accept parameters " + p + ", " + N;
			return -1;
		}
*/
	}
	
	
	private boolean alreadyContainsSite(int r, int numSites) {
		for (int i = 0; i < numSites; i++)
			if (randomSite[i] == r)
				return true;
		return false;
	}
	
	
	private void selectRandomSitesInRange(int site, int numSites) {
		assert (numSites < 12);
		
		for (int i = 0; i < numSites; i++) {
			int r;
			do {
				r = (int) (Math.random()*R + 1 + site);
				r = (r+L) % L;
			} while (r == site || spins.get(r) != direction || alreadyContainsSite(r, i));
			randomSite[i] = r;
		}
	}
	
	
	private void sanityCheck() {
		int sum = 0;
		for (int size = 0; size <= L; size++) {
			sum += size*numClusters[size];
		}
		for (int site = 0; site < L; site++) {
			if (spins.get(site) == direction)
				sum--;
		}
		assert(sum == 0);
	}
}

