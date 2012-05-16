package kip.javasim.molecular;


public class Atom3d<Tag> extends Pt3d {
	final Tag tag;
	final int wrapX, wrapY, wrapZ;
	
	public Atom3d(double x, double y, double z) {
		super(x,y,z);
		wrapX = wrapY = wrapZ = 0;
		tag = null;
	}
}
