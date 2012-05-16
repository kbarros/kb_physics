package kip.javasim.md.apps.shaker;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;

import kip.javasim.md.Particle;
import kip.javasim.md.ParticleContext;


class SimulationTrajectory implements AbstractTrajectory {
	private ArrayList<Snapshot> snapshots;
	private double t_i, t_f;
	
	public SimulationTrajectory(String dir) {
		String[] fs = (new File(dir)).list();
		if (fs == null)
			System.err.println("Directory '" + dir + "' doesn't exist.");
		
		snapshots = new ArrayList<Snapshot>();
		for (String f : fs) {
			if (!f.substring(0, 2).equals("t="))
				continue;
			double time = Double.valueOf(f.substring(2));
			snapshots.add(new Snapshot(dir+File.separator+f, time));
		}
	    Collections.sort(snapshots);
	    
	    t_i = snapshots.get(0).time;
	    t_f = snapshots.get(snapshots.size()-1).time;
	}
	
	public double startTime() {
		return t_i;
	}
	
	public double endTime() {
		return t_f;
	}
	
	public Particle[] get(double t) {
		// perform a binary search for snapshot at nearest time
		int lb = 0; // lower bound
		int ub = snapshots.size() - 1; // upper bound
		while (ub - lb > 1) {
			int i = (lb + ub) / 2;
			if (time(i) > t)
				ub = i;
			else
				lb = i;
		}
		double dt1 = t - time(lb);
		double dt2 = time(ub) - t;
		int i = (dt1 < dt2) ? lb : ub;
		return snapshots.get(i).particles(); 
	}
	
	public ParticleContext getContext() {
		return snapshots.get(0).particles()[0].tag.pc;
	}
	
	private double time(int i) {
		return snapshots.get(i).time();
	}
	
	class Snapshot implements Comparable<Snapshot> {
		private double time;
		private String fname;
		private Particle[] ps;
		
		public Snapshot(String fname, double time) {
			this.fname = fname;
			this.time = time;
		}
		
		public double time() {
			return time;
		}
		
		public Particle[] particles() {
			if (ps == null)
				ps = ParticleContext.readParticles(fname);
			return ps;
		}
		
		public int compareTo(Snapshot that) {
			return Double.compare(time, that.time);
		}
	}
}
