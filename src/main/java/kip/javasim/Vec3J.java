package kip.javasim;

import static java.lang.Math.*;

// TODO replace with scikit.vecmath.Vector3d
public class Vec3J implements Cloneable {
	public double x, y, z;
	public Vec3J(double x, double y, double z) {
		this.x = x;
		this.y = y;
		this.z = z;
	}
	
	public Vec3J clone(Vec3J that) {
		try {
			return (Vec3J)super.clone();
		} catch (Exception e) {
			return null;
		}
	}
	
	public double dot(Vec3J that) {
		return x*that.x + y*that.y + z*that.z;
	}
	
	public double norm2() {
		return this.dot(this);
	}
	
	public double norm() {
		return Math.sqrt(norm2());
	}
	
	public Vec3J normalize() {
		if (norm2() == 0)
			throw new IllegalArgumentException("Can't normalize zero vector.");
		else
			return scale(1/norm());
	}
	
	public Vec3J projectOnto(Vec3J that) {
		// v1*v2 v2 / v2^2
		return that.scale(this.dot(that) / that.norm2());
	}
	
	public Vec3J scale(double a) {
		return new Vec3J(x*a, y*a, z*a);
	}
	
	public Vec3J plus(Vec3J that) {
		return new Vec3J(x+that.x, y+that.y, z+that.z);
	}
	
	public Vec3J minus(Vec3J that) {
		return new Vec3J(x-that.x, y-that.y, z-that.z);		
	}
	
	public Vec3J cross(Vec3J v) {
		return new Vec3J(y*v.z - z*v.y, -(x*v.z - z*v.x), x*v.y - y*v.x);
	}
	
	public Vec3J rotateZ(double a) {
		double c = cos(a);
		double s = sin(a);
		return new Vec3J(c*x-s*y, s*x+c*y, z);
	}
	
	public String toString() {
		return "("+x+","+y+","+z+")";
	}
}
