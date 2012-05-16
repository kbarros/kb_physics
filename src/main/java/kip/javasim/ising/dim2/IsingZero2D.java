package kip.javasim.ising.dim2;

public class IsingZero2D extends Ising2D {	
	int NOT_ACTIVE = -1;
	int activityIndex[]; // either -1 (inactive) or index into activeSites[]
	int activeSites[]; // active site indices
	int activeSitesCnt; // number of active sites
	
	public IsingZero2D(int seed, int _L1, int _L2, double _T, boolean _openBoundary) {
		super(seed, _L1, _L2, _T, _openBoundary);
		if (T != 0)
			throw new IllegalArgumentException("Zero temperature only.");
		activityIndex = new int[N];
		activeSites = new int[N];
		activeSitesCnt = 0;
		for (int i = 0; i < N; i++) {
			activityIndex[i] = NOT_ACTIVE;
			updateSiteActivation(i);
		}
	}
	
	private void updateSiteActivation(int i) {
		boolean shouldBeActive = (spin[i] * neighborSum(i)) <= 0;
		
		if (activityIndex[i] == NOT_ACTIVE) {
			// add site to active list if necessary
			if (shouldBeActive) {
				activityIndex[i] = activeSitesCnt;
				activeSites[activeSitesCnt] = i;
				activeSitesCnt++;
			}
		}
		else {
			// remove site from active list if necessary
			if (!shouldBeActive) {
				activeSitesCnt--;
				int availableIdx = activityIndex[i];
				int siteToMove = activeSites[activeSitesCnt];
				activityIndex[siteToMove] = availableIdx;
				activityIndex[i] = NOT_ACTIVE;
				activeSites[availableIdx] = siteToMove;
			}
		}
	}
	
	private void singleStep() {
		if (activeSitesCnt == 0)
			return;
		int i = activeSites[random.nextInt(activeSitesCnt)];
		double dE = 2*spin[i]*neighborSum(i);
		if (dE <= 0) {
			spin[i] = -spin[i];
			updateSiteActivation(i);
			updateSiteActivation(i % L1 == 0    ? i+(L1-1) : i-1); // left
			updateSiteActivation(i % L1 == L1-1 ? i-(L1-1) : i+1); // right
			updateSiteActivation((i-L1+N) % N); // down
			updateSiteActivation((i+L1  ) % N); // up
		}
	}

	
	public void step(double mcs) {
		double t = 0;
		while (t < mcs && activeSitesCnt > 0) {
			t += 1. / activeSitesCnt;  
			singleStep();
		}
		if (activeSitesCnt == 0)
			time += mcs;
		else
			time += t;
	}

}
