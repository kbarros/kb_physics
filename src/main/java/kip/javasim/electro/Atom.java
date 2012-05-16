package kip.javasim.electro;

import scikit.util.Utilities;


public class Atom {
	final Maggs maggs;
	final int tag;
	final double charge;
	final double x, y, z;
	final int wrapX, wrapY, wrapZ;
	
	public Atom(Maggs maggs, int tag, double charge, double x, double y, double z, int wrapX, int wrapY, int wrapZ) {
		this.maggs = maggs;
		this.tag = tag;
		this.charge = charge;
		this.x = x;
		this.y = y;
		this.z = z;
		this.wrapX = wrapX;
		this.wrapY = wrapY;
		this.wrapZ = wrapZ;
		
		double L = maggs.L;
		while (x < 0)  {wrapX--; x += L;};
		while (x >= L) {wrapX++; x -= L;};
		while (y < 0)  {wrapY--; y += L;};
		while (y >= L) {wrapY++; y -= L;};
		while (z < 0)  {wrapZ--; z += L;};
		while (z >= L) {wrapZ++; z -= L;};	
	}
	
	public Atom(Maggs maggs, int tag, double charge) {
		this(maggs, tag, charge, 0., 0., 0., 0, 0, 0);
	}
	
	public Atom displace(double dx, double dy, double dz) {
		return new Atom(maggs, tag, charge, x+dx, y+dy, z+dz, wrapX, wrapY, wrapZ);
	}
		
	public double shortestDistance(Atom that) {
		double dx = Utilities.periodicOffset(maggs.L, that.x - x);
		double dy = Utilities.periodicOffset(maggs.L, that.y - y);
		double dz = Utilities.periodicOffset(maggs.L, that.z - z);
		
		return Math.sqrt(dx*dx + dy*dy + dz*dz);
	}
}
