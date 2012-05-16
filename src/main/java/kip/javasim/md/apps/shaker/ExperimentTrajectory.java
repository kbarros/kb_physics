package kip.javasim.md.apps.shaker;

import static scikit.util.FileUtil.readFloatLittleEndian;
import static scikit.util.FileUtil.readIntLittleEndian;

import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;

import kip.javasim.md.Particle;
import kip.javasim.md.ParticleContext;

public class ExperimentTrajectory implements AbstractTrajectory {
	double t_i, t_f;
	public Particle[][] snapshots;
	
	@SuppressWarnings("unused")
  public ExperimentTrajectory(String filename) {
		try {
		    DataInputStream dis = new DataInputStream(new FileInputStream(filename));
		    
		    readIntLittleEndian(dis); // unknown meaning
		    int dim = readIntLittleEndian(dis); // spatial dimensions
		    int n = readIntLittleEndian(dis); // columns of data
		    int m = readIntLittleEndian(dis); // rows of data
		    int prec = readIntLittleEndian(dis); // precision (4 for float, 8 for double)
		    int size = readIntLittleEndian(dis); // m*n
		    for (int i = 0; i < 6; i++)
		    	readIntLittleEndian(dis); // these should be zero
		    assert(dim == 2);
		    assert(prec == 4);
		    assert(size == m*n);
		    
		    System.out.println("m " + m);
		    for (int i = 0; i < m; i++) {
		    	double x = readFloatLittleEndian(dis);
		    	double y = readFloatLittleEndian(dis);
		    	double brightness = readFloatLittleEndian(dis);
		    	double radius = readFloatLittleEndian(dis);
		    	double time = readFloatLittleEndian(dis);
		    	double id = readFloatLittleEndian(dis);
		    	x += y+brightness+radius+time+id;
//		    	System.out.println(time + " "+ id);
		    }
		    
		    // read context
//		    double L = dis.readDouble();
//		    Type type = dis.readBoolean() ? Type.Torus2D : Type.Disk2D;
//		    ParticleContext pc = new ParticleContext(L, type);
		    
		    // read tags
//		    Map<Integer, ParticleTag> tags = new TreeMap<Integer, ParticleTag>();
//		    int ntags = dis.readInt();
		    /*
		    for (int i = 0; i < ntags; i++) {
		    	ParticleTag tag = new ParticleTag(dis.readInt());
		    	tag.pc = pc;
		    	tag.mass = dis.readDouble();
		    	tag.radius = dis.readDouble();
		    	tag.color = new Color(dis.readInt(), dis.readInt(), dis.readInt(), 128);
		    	tags.put(tag.id, tag);
		    }
		    
		    // read particles
			ret = new Particle[dis.readInt()];
			for (int i = 0; i < ret.length; i++) {
				ret[i] = new Particle();
			    ret[i].tag = tags.get(dis.readInt());
			    ret[i].x = dis.readDouble();
			    ret[i].y = dis.readDouble();
			    ret[i].z = dis.readDouble();
			}*/
		    
		    dis.close();
		}
		catch (IOException e) {
			System.err.println ("Error reading from '" + filename + "'");
		}
	}
	
	public double startTime() {
		return 0;
	}
	public double endTime() {
		return 0;
	}
	public Particle[] get(double t) {
		return null;
	}
	public ParticleContext getContext() {
		return null;
	}
}
