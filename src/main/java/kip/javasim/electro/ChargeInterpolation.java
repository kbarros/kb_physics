package kip.javasim.electro;

import scikit.util.Pair;


abstract public class ChargeInterpolation {
	static class Ret {
		int[] sites;
		double[] densityChange;
		int[] links;
		double[] electricChange;
		
		public Ret(int[] sites, double[] densityChange, int[] links, double[] electricChange) {
			this.sites = sites;
			this.densityChange = densityChange;
			this.links = links;
			this.electricChange = electricChange;
		}
	}
	
	final Maggs maggs;
	
	ChargeInterpolation(Maggs maggs) {
		this.maggs = maggs;
	}
	
	// Interpolates the charge of an atom onto the lattice.  Return the relevant sites and their charge
	// density
	abstract public Pair<int[], double[]> interpolateCharge(Atom a);
	
	// Calculates updates to the interpolated charge and electric field of an atom that is displaced.
	// Returns the affected site and links indices and their deltas. 
	public Ret siteAndLinkDeltas(Atom a, double dx, double dy, double dz) {
		if (dx >= 0 && dy == 0 && dz == 0)
			return deltasX(a.charge, dx, a.x, a.y, a.z);
		else if (dx < 0 && dy == 0 && dz == 0)
			return deltasX(-a.charge, -dx, a.x+dx, a.y, a.z);
		else if (dx == 0 && dy >= 0 && dz == 0)
			return deltasY(a.charge, dy, a.x, a.y, a.z);
		else if (dx == 0 && dy < 0 && dz == 0)
			return deltasY(-a.charge, -dy, a.x, a.y+dy, a.z);
		else if (dx == 0 && dy == 0 && dz >= 0)
			return deltasZ(a.charge, dz, a.x, a.y, a.z);
		else if (dx == 0 && dy == 0 && dz < 0)
			return deltasZ(-a.charge, -dz, a.x, a.y, a.z+dz);
		else {
			System.err.println("Steps must be along a single axis");
			return null;
		}
	}
	
	abstract Ret deltasX(double charge, double dx, double x, double y, double z);
	abstract Ret deltasY(double charge, double dy, double x, double y, double z);
	abstract Ret deltasZ(double charge, double dz, double x, double y, double z);
	
	
	
	static public class Lattice extends ChargeInterpolation {
		public Lattice(Maggs maggs) {
			super(maggs);
		}
		
		public Pair<int[], double[]> interpolateCharge(Atom a) {
			int i, j, k;
			i = maggs.getNearestCoordIndex(a.x);
			j = maggs.getNearestCoordIndex(a.y);
			k = maggs.getNearestCoordIndex(a.z);
			int site = maggs.getSiteAtIndices(i, j, k);
			return new Pair<int[],double[]>(new int[]{site}, new double[]{a.charge});
		}
		
		Ret deltasX(double charge, double dx, double x, double y, double z) {
			if (dx > 1)
				throw new IllegalArgumentException("Step size too big");
			int i1 = maggs.getNearestCoordIndex(x);
			int i2 = maggs.getNearestCoordIndex(x+dx);
			int j = maggs.getNearestCoordIndex(y);
			int k = maggs.getNearestCoordIndex(z);
			int site1 = maggs.getSiteAtIndices(i1, j, k);
			int site2 = maggs.getSiteAtIndices(i2, j, k);
			return deltasHelper(charge, site1, site2, 3*site1+maggs.XDir);
		}
		
		Ret deltasY(double charge, double dy, double x, double y, double z) {
			if (dy > 1)
				throw new IllegalArgumentException("Step size too big");
			int i = maggs.getNearestCoordIndex(x);
			int j1 = maggs.getNearestCoordIndex(y);
			int j2 = maggs.getNearestCoordIndex(y+dy);
			int k = maggs.getNearestCoordIndex(z);
			int site1 = maggs.getSiteAtIndices(i, j1, k);
			int site2 = maggs.getSiteAtIndices(i, j2, k);
			return deltasHelper(charge, site1, site2, 3*site1+maggs.YDir);
		}
		
		Ret deltasZ(double charge, double dz, double x, double y, double z) {
			if (dz > 1)
				throw new IllegalArgumentException("Step size too big");
			int i = maggs.getNearestCoordIndex(x);
			int j = maggs.getNearestCoordIndex(y);
			int k1 = maggs.getNearestCoordIndex(z);
			int k2 = maggs.getNearestCoordIndex(z+dz);
			int site1 = maggs.getSiteAtIndices(i, j, k1);
			int site2 = maggs.getSiteAtIndices(i, j, k2);
			return deltasHelper(charge, site1, site2, 3*site1+maggs.ZDir);
		}
		
		private Ret deltasHelper(double charge, int site1, int site2, int link) {
			if (site1 == site2) {
				System.err.println("Ineffective move!");
				return new Ret(
						new int[]{},
						new double[]{},
						new int[]{},
						new double[]{});
			}
			else {
				return new Ret(
						new int[]{site1, site2},
						new double[]{-charge, charge},
						new int[]{link},
						new double[]{- charge / maggs.dielectric});
			}

		}
	}
}
