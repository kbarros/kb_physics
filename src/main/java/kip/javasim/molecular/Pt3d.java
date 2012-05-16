package kip.javasim.molecular;


public class Pt3d {
	final double x, y, z;
	
	public Pt3d(double x, double y, double z) {
		this.x = x;
		this.y = y;
		this.z = z;
	}
	
	public Pt3d displace(double dx, double dy, double dz) {
		return new Pt3d(x+dx, y+dy, z+dz);
	}
}

/*
public class Atom {
	final Maggs maggs;
	final int tag;
	final double charge;
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
	
*/