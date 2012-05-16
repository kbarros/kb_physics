package kip.javasim.electro;

import java.util.Random;

import scikit.util.Pair;


public class Maggs {
	// parameters
	double L = 5;
	int cols = 5;
	double dx = L/cols;
	double dx3 = dx*dx*dx;
	double dielectric = 1.0;
	double Emag = 1.0;
	double beta = 1;
	double mu = 0.75/dx; // inverse length scale of Yukawa interactions 
	
	Random rand;
	HashedArrayList<Atom> atoms;
	BinnedAtoms<Atom> bins;
	ChargeInterpolation interp;
	
	// mcs acceptance statistics
	int rotatePlaquetteAccept = 0;
	int rotatePlaquetteTrials = 0;
	int shiftNetEAccept = 0;
	int shiftNetETrials = 0;
	int moveChargeAccept = 0;
	int moveChargeTrials = 0;
	
	// E field layout
	int _numSites = cols*cols*cols; 
	int _numLinks = _numSites*3;
	int XDir = 0;
	int YDir = 1;
	int ZDir = 2;
	
	// layout of E field is (from slowest to fastest index)
	// [z] [y] [x] [link_idx] 
	// where 'link_idx' chooses from the set of three outgoing links {XDir YDir ZDir}
	// in the "positive" direction
	double E[] = new double[_numLinks];
	
	// psi field generates Yukawa interactions, as described in Rottler/Maggs JCP 2004
	double psi[] = new double[_numSites];
	
	// charge at each site (used for testing only)
	double rho[] = new double[_numSites];
	
	// energy in E and psi fields
	double fieldEnergy;
	
	
	public Maggs() {
		 rand = new Random(0);
		 atoms = new HashedArrayList<Atom>();
		 bins = new BinnedAtoms<Atom>(L, cols);
		 interp = new ChargeInterpolation.Lattice(this);
		 
		 // system is initially empty
		 for (int link = 0; link < _numLinks; link++) {
			 E[link] = 0;
		 }
		 for (int site = 0; site < _numSites; site++) {
			 psi[site] = 0;
			 rho[site] = 0;
		 }
	}
	
	public void createAtoms(Atom... newAtoms) {
		// check charge neutrality of collection of atoms
		double netCharge = 0;
		for (Atom a : newAtoms)
			netCharge += a.charge;
		if (Math.abs(netCharge) > 1e-8)
			throw new IllegalArgumentException("Attempt to create unbalanced charge distribution");
		
		// create each atom by translating from (0,0,0) to desired locations.
		// Gauss' law is satisfied by construction
		for (Atom a : newAtoms) {
			// move atom to (0,0,0)
			Atom ap = a.displace(-a.x, -a.y, -a.z);			
			// add atom to data structures
			atoms.add(ap);
			bins.addAtom(ap);
			// move atom to desired location in numSteps
			int numSteps = cols;
			double dx = a.x/numSteps, dy = a.y/numSteps, dz = a.z/numSteps;
			for (int i = 0; i < numSteps; i++) {
				ap = displaceAtomAndUpdateFields(ap, dx, 0, 0);
				ap = displaceAtomAndUpdateFields(ap, 0, dy, 0);
				ap = displaceAtomAndUpdateFields(ap, 0, 0, dz);
			}
		}
	}
	
	public Atom atom(int i) {
		return atoms.get(i);
	}
	
	public void doMCSSweep() {
		for (int cnt = 0; cnt < _numLinks; cnt++) {
			int i = rand.nextInt(cols);
			int j = rand.nextInt(cols);
			int k = rand.nextInt(cols);
			rotatePlaquette(i, j, k);
			// moveCharge(i, j, k);
		}
		shiftNetE();
	}
	
	int getNearestCoordIndex(double x) {
		return ((int)Math.round(x/dx)) % cols;
	}
	
	int getIndexXAtSite(int site) {
		return site % cols;
	}
	
	int getIndexYAtSite(int site) {
		return (site/cols) % cols;
	}
	
	int getIndexZAtSite(int site) {
		return site / (cols*cols);
	}
	
	int getSiteAtIndices(int i, int j, int k) {
		i = (i + cols) % cols;
		j = (j + cols) % cols;
		k = (k + cols) % cols;
		return ((k)*cols + j)*cols + i;
	}
	
	int getLinkAtIndices(int i, int j, int k, int dir) {
		i = (i + cols) % cols;
		j = (j + cols) % cols;
		k = (k + cols) % cols;
		return (((k)*cols + j)*cols + i)*3 + dir;
	}
	
	Atom displaceAtomAndUpdateFields(Atom a, double dx, double dy, double dz) {
		ChargeInterpolation.Ret ret = interp.siteAndLinkDeltas(a, dx, dy, dz);
		for (int i = 0; i < ret.links.length; i++)
			E[ret.links[i]] += ret.electricChange[i];
		for (int i = 0; i < ret.sites.length; i++)
			rho[ret.sites[i]] += ret.densityChange[i];
		Atom ap = a.displace(dx, dy, dz);
		atoms.replace(a, ap);
		bins.removeAtom(a);
		bins.addAtom(ap);
		return ap;
	}
	
	double getChargeAtSite(int site) {
		int i = getIndexXAtSite(site);
		int j = getIndexYAtSite(site);
		int k = getIndexZAtSite(site);
		double sum =
			E[getLinkAtIndices(i, j, k, XDir)] - E[getLinkAtIndices(i-1, j, k, XDir)] +
			E[getLinkAtIndices(i, j, k, YDir)] - E[getLinkAtIndices(i, j-1, k, YDir)] +
			E[getLinkAtIndices(i, j, k, ZDir)] - E[getLinkAtIndices(i, j, k-1, ZDir)];
		return sum*dielectric;
	}
	
	
	void consistencyCheck() {
		// Make sure that the charge is self consistent
		double[] rhoTest = new double[_numSites];
		for (Atom a : atoms) {
			Pair<int[],double[]> charges = interp.interpolateCharge(a);
			for (int i = 0; i < charges.fst().length; i++) {
				int site = charges.fst()[i];
				double q = charges.snd()[i]; 
				rhoTest[site] += q;
			}
		}
		for (int site = 0; site < _numSites; site++) {
			double q1 = rhoTest[site];
			double q2 = rho[site];
			double q3 = getChargeAtSite(site);
			if (Math.abs(q1 - q2) > 1e-8 || Math.abs(q2 - q3) > 1e-8) {
				System.err.println("Charges (" + q1+" "+q2+" "+q3+") at site "+site);
			}
		}
		
		// Make sure that the energy is self consistent
		double fieldEnergy2 = calculateFieldEnergy();
		if (Math.abs(fieldEnergy - fieldEnergy2) > 1e-8) {
			System.err.println("Invalid energy ("+fieldEnergy+" "+fieldEnergy2);
		}
	}
	
	void checkAcceptanceRatios() {
		double r1 = rotatePlaquetteAccept/(double)rotatePlaquetteTrials;
		double r2 = shiftNetEAccept/(double)shiftNetETrials;
		double r3 = moveChargeAccept/(double)moveChargeTrials;
		System.out.println("ratios: " + r1 + " " + r2 + " " + r3);
	}
	
	// Returns the energy contained in the electric field E and scalar field psi 
	double calculateFieldEnergy() {
		double U = 0;
		for (int link = 0; link < _numLinks; link++) {
			// energy due to x/y/z component of electric field
			U += 0.5*dielectric*dx3*E[link]*E[link];
		}
		for (int site = 0; site < _numSites; site++) {
			// (\eps_0 / 2) (grad \psi)^2
			int i = getIndexXAtSite(site);
			int j = getIndexYAtSite(site);
			int k = getIndexZAtSite(site);
			double dpx = (psi[getSiteAtIndices(i+1, j, k)] - psi[site]) / dx;
			double dpy = (psi[getSiteAtIndices(i, j+1, k)] - psi[site]) / dx;
			double dpz = (psi[getSiteAtIndices(i, j, k+1)] - psi[site]) / dx;
			double grad_psi_2 = dpx*dpx + dpy*dpy + dpz*dpz;
			U += 0.5 * dielectric * dx3 * grad_psi_2;
			
			// (\eps_0 / 2) \mu^2 \psi^2
			U += 0.5 * dielectric * dx3 * mu*mu * psi[site]*psi[site]; 
			
			// - \rho \psi
			U += - dx3 * psi[site] * rho[site];
		}
		return U;
	}
	
	
	// Returns the energy cost to add 'dE' to the electric field at 'link'
	double deltaElectricEnergy(int link, double dE) {
		return 0.5 * dielectric * dx3 * (2*E[link]*dE + dE*dE);
	}
	
	// Returns the energy cost to add 'dRho' to the charge at 'site'.
	double deltaChargeEnergy(int site, double dRho) {
		return - dx3 * psi[site] * dRho;
	}
	
	// Returns the energy cost to add 'dPsi' to the scalar field at 'site'
	double deltaPsiEnergy(int site, double dPsi) {
		int i = getIndexXAtSite(site);
		int j = getIndexYAtSite(site);
		int k = getIndexZAtSite(site);
		
		double ret = 0;

		// (\eps_0 / 2) (grad \psi)^2
		//
		// note: relevant contribution to [ 1/2 (grad \psi)^2 ] is, in each direction,
		//   1/2 [ ((psi1-psi0)/dx)^2 + ((psi2-psi1)/dx)^2 ]
		//     = ( psi1^2 - psi0 psi1 - ps1 ps2 ) / dx^2
		double xpart = - psi[getSiteAtIndices(i+1, j, k)] - psi[getSiteAtIndices(i-1, j, k)];
		double ypart = - psi[getSiteAtIndices(i, j+1, k)] - psi[getSiteAtIndices(i, j-1, k)];
		double zpart = - psi[getSiteAtIndices(i, j, k+1)] - psi[getSiteAtIndices(i, j, k-1)];
		ret += dielectric * dx3 * (3*(2*psi[site]*dPsi + dPsi*dPsi) + xpart + ypart + zpart) / (dx * dx);
		
		// (\eps_0 / 2) \mu^2 \psi^2
		ret += 0.5 * dielectric * dx3 * mu*mu * (2*psi[site]*dPsi + dPsi*dPsi);

		// - \rho \psi
		ret += - dx3 * getChargeAtSite(site) * dPsi;
		
		return ret;
	}
	
	// Returns the energy cost to displace a particle according to 'disp' 
	double displacementDeltaEnergy(ChargeInterpolation.Ret disp) {
		double ret = 0;
		for (int i = 0; i < disp.links.length; i++)
			ret += deltaElectricEnergy(disp.links[i], disp.electricChange[i]);
		for (int i = 0; i < disp.sites.length; i++)
			ret += deltaChargeEnergy(disp.sites[i], disp.densityChange[i]);
		return ret;
	}

	void rotatePlaquette(int x, int y, int z) {
		
		//
		// Links on a plaquette have the structure
		//
		//      2
		//    ----->
		//   ^      ^
		// 1 |      | 3
		//   |      |
		//  (X)---->
		//      0
		//
		// where the site (x,y,z) is located at position (X)
		//
		
		// randomly select from 1 of 3 possible plaquettes
		int i0, i1, i2, i3;
		switch (rand.nextInt(3)) {
		// x-y plaquette
		case 0:
			i0 = getLinkAtIndices(x,y,z,XDir);
			i1 = getLinkAtIndices(x,y,z,YDir);
			i2 = getLinkAtIndices(x,y+1,z,XDir);
			i3 = getLinkAtIndices(x+1,y,z,YDir);
			break;
		// x-z plaquette
		case 1:
			i0 = getLinkAtIndices(x,y,z,XDir);
			i1 = getLinkAtIndices(x,y,z,ZDir);
			i2 = getLinkAtIndices(x,y,z+1,XDir);
			i3 = getLinkAtIndices(x+1,y,z,ZDir);
			break;
		// y-z plaquette
		case 2:
			i0 = getLinkAtIndices(x,y,z,YDir);
			i1 = getLinkAtIndices(x,y,z,ZDir);
			i2 = getLinkAtIndices(x,y,z+1,YDir);
			i3 = getLinkAtIndices(x,y+1,z,ZDir);
			break;
		default:
			throw new IllegalStateException();
		}
		
		// generate trial change to electric field
		double dE = Emag*(rand.nextDouble()-0.5);
		double dE0 = +dE;
		double dE1 = -dE;
		double dE2 = -dE;
		double dE3 = +dE;
		
		// calculate change in energy for the plaquette rotation
		double deltaEnergy = deltaElectricEnergy(i0,dE0) + deltaElectricEnergy(i1,dE1) + deltaElectricEnergy(i2,dE2) + deltaElectricEnergy(i3,dE3);
		
		// accept new state according to the metropolis criterion
		if (deltaEnergy < 0 || rand.nextDouble() < Math.exp(-beta*deltaEnergy)) {
			E[i0] += dE0;
			E[i1] += dE1;
			E[i2] += dE2;
			E[i3] += dE3;
			fieldEnergy += deltaEnergy;
			rotatePlaquetteAccept++;
		}
		rotatePlaquetteTrials++;
	}
	
	void shiftNetE() {
		double norm = Emag/Math.sqrt(_numLinks);
		double dEx = norm*(rand.nextDouble()-0.5);
		double dEy = norm*(rand.nextDouble()-0.5);
		double dEz = norm*(rand.nextDouble()-0.5);
		
		double deltaEnergy = 0;
		for (int i = 0; i < _numSites; i++) {
			deltaEnergy += deltaElectricEnergy(3*i+XDir, dEx);
			deltaEnergy += deltaElectricEnergy(3*i+YDir, dEy);
			deltaEnergy += deltaElectricEnergy(3*i+ZDir, dEz);
		}
		
		if (deltaEnergy < 0 || rand.nextDouble() < Math.exp(-beta*deltaEnergy)) {
			for (int i = 0; i < _numSites; i++) {
				E[3*i+0] += dEx;
				E[3*i+1] += dEy;
				E[3*i+2] += dEz;
			}
			fieldEnergy += deltaEnergy;
			shiftNetEAccept++;
		}
		shiftNetETrials++;
	}
	
	
	void moveAtom(Atom a) {
		
	}
	
//	void moveCharge(int x, int y, int z) {
//		int charge1, charge2;
//		int i; // link from charge1 to charge2
//		
//		switch (rand.nextInt(6)) {
//		case 0:
//			charge1 = getChargeAtSite(x-1,y,z);
//			charge2 = getChargeAtSite(x, y, z);
//			i = getLinkAtIndices(x-1,y,z,XDir);
//			break;
//		case 1:
//			charge1 = getChargeAtSite(x, y, z);
//			charge2 = getChargeAtSite(x+1,y,z);
//			i = getLinkAtIndices(x,y,z,XDir);
//			break;
//		case 2:
//			charge1 = getChargeAtSite(x,y-1,z);
//			charge2 = getChargeAtSite(x, y, z);
//			i = getLinkAtIndices(x,y-1,z,YDir);
//			break;
//		case 3:
//			charge1 = getChargeAtSite(x, y, z);
//			charge2 = getChargeAtSite(x,y+1,z);
//			i = getLinkAtIndices(x,y,z,YDir);
//			break;
//		case 4:
//			charge1 = getChargeAtSite(x,y,z-1);
//			charge2 = getChargeAtSite(x, y, z);
//			i = getLinkAtIndices(x,y,z-1,ZDir);
//			break;
//		case 5:
//			charge1 = getChargeAtSite(x, y, z);
//			charge2 = getChargeAtSite(x,y,z+1);
//			i = getLinkAtIndices(x,y,z,ZDir);
//			break;
//		default: throw new IllegalStateException();
//		}
//		
//		
//		// dE is the trial change to E on the link from charge1 to charge2.
//		// positive dE corresponds to _increasing_ charge1 and _decreasing_ charge2
//		double dE;
//		if (charge1 != 0 && charge2 == 0)
//			dE = - charge1 / dielectric;
//		else if (charge1 == 0 && charge2 != 0)
//			dE = + charge2 / dielectric;
//		else
//			return;
//		
//		double deltaEnergy = linkDeltaEnergy(i,dE);
//		
//		// accept new state according to the metropolis criterion
//		if (deltaEnergy < 0 || rand.nextDouble() < Math.exp(-beta*deltaEnergy)) {
//			E[i] += dE;
//			moveChargeAccept++;
//		}
//		moveChargeTrials++;
//	}
}
