package kip.javasim.electro.apps;

import kip.javasim.electro.Atom;
import kip.javasim.electro.Maggs;
import scikit.dataset.Histogram;

public class PairPotentialApp {
	public static void main(String[] args) {
		Maggs maggs = new Maggs();
		
		Atom a1 = new Atom(maggs, 0, +1.).displace(0,0,0);
		Atom a2 = new Atom(maggs, 0, -1.).displace(0,0,0);
		maggs.createAtoms(a1, a2);
		
		Histogram pairDistance = new Histogram(0.2);
		pairDistance.setNormalizing(true);
		
		for (int i = 0; i < 100; i++)
			for (int j = 0; j < 10; j++) {
				maggs.doMCSSweep();
			pairDistance.accum(getPairDistance(maggs));
			
			// checkE();
			// checkAcceptanceRatios();
		}
		
		scikit.util.FileUtil.dumpColumns("energy", pairDistance.copyData().columns());
	}
	
	
	public static double getPairDistance(Maggs maggs) {
		Atom a1 = maggs.atom(0);
		Atom a2 = maggs.atom(1);
		return a1.shortestDistance(a2);
	}
}
