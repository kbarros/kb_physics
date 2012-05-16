package kip.javasim.ising;


abstract public class RewindableDynamics implements Cloneable {
	private RewindableDynamics old;
	private double memoryTime = Double.MAX_VALUE;
	protected double time, dt;
	
	
	abstract protected void _step();
	
	
	public RewindableDynamics clone() {
		try {
			return (RewindableDynamics)super.clone();
		} catch (CloneNotSupportedException e) {
			return null;
		}
	}
	
	
	public void setMemoryTime(double memoryTime) {
		this.memoryTime = memoryTime;		
	}
	
	
	public RewindableDynamics simulationAtTime(double t) {
		if (t < time()) {
			return (old == null) ? null : old.simulationAtTime(t);
		}
		else {
			assert(time() <= t);
			RewindableDynamics c = clone();
			c.runUntil(t);
			return c;
		}
	}
	
	
	public void step() {
		if (old != null) {
			assert (time() >= old.time());
			if (old.old != null) {
				assert (old.time() >= old.old.time());
				assert (old.old.old == null);
			}
		}
		if (old == null)
			old = clone();
		_step();
		time += dt; // BUG: some error here
		if (time() - old.time() > memoryTime) {
			old = clone();
			old.old.old = null; // cut off history to avoid memory leaks
		}
	}
	
	
	public void runUntil(double t) {
		while (time() < t)
			step();
	}
	
	
	public double time() {
		return time;
	}
	
	
	public void resetTime() {
		time = 0;
		old = null;
	}
}
