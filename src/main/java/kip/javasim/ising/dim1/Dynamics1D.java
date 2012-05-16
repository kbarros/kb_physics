package kip.javasim.ising.dim1;

import scikit.jobs.Job;
import scikit.jobs.params.Parameters;
import kip.javasim.Random;
import static java.lang.Math.*;


abstract public class Dynamics1D implements Cloneable {
	private Dynamics1D		old;
	private double			memoryTime;
	
	public Random			random = new Random();
	public int N, R, dx;
	public double time, dt, h;
	
	static final double NUCLEATION_CUTOFF = 0;
    static final double NUCLEATION_TIME_ACCURACY = 0.05;
    static final double NUCLEATION_WAIT_TIME = 40;
    
    
	public Dynamics1D clone() {
        try {
            Dynamics1D c = (Dynamics1D)super.clone();
            c.random = (Random)random.clone();
            return c;
        } catch (Exception e) {
            return null;
        }
	}
	
	
	public Dynamics1D simulationAtTime(double t) {
		if (t < time()) {
			return (old == null) ? null : old.simulationAtTime(t);
		}
		else {
			assert(time() <= t);
			Dynamics1D c = clone();
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
		Job.yield();
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
	
	
	public void initialize(Parameters params) {
		time = 0;
		old = null;
		
		memoryTime = params.fget("Memory time", Double.POSITIVE_INFINITY);
		random.setSeed(params.iget("Random seed", 0));
		
		setParameters(params);
	}
	
	
	public void setParameters(Parameters params) {
		dt = params.fget("dt");
		h  = params.fget("h", 0);		
	}
	
	
	public double[] copyField() {
		double ret[] = new double[N/dx];
		for (int i = 0; i < N/dx; i++)
			ret[i] = fieldElement(i);
		return ret;
	}
	
	
	abstract public double magnetization();
	abstract public void randomizeField(double m);
	abstract public void setField(double m);
	abstract public double fieldElement(int i);
	abstract protected void _step(); // step without saving "old" sim copies
    
    
    // -- NUCLEATION ----------------------------------------------------------------
    //
    
    public scikit.dataset.Function saddleProfile() {
        return null;
    }
    
	public boolean nucleated() {
		for (int i = 0; i < N/dx; i++)
			if (fieldElement(i) > NUCLEATION_CUTOFF)
				return true;
		return false;
	}
    
    public static int maxIndex(double[] a) {
        int imax = 0;
        for (int i = 0; i < a.length; i++)
			imax = a[i] > a[imax] ? i : imax;
        return imax;
    }
    
    // return time and location
	public double[] nucleationTimeAndLocation(double overshootEstimate) {
        double[] a = null;
        double loc = dx*maxIndex(copyField());
        
        double t_hi = time();
        double t_lo = max(time()-2*overshootEstimate, 0);
        double t_md = (t_hi + t_lo) / 2;
        
        // binary search for nucleation time
        while (t_hi - t_lo > NUCLEATION_TIME_ACCURACY) {
            a = simulationAtTime(t_md).nucleationProbabilityAndLocation(loc);
            if (a[0] < 0.5) {
                t_lo = t_md;
            }
            else {
                t_hi = t_md;
                loc = a[1];
            }
            t_md = (t_hi + t_lo) / 2;
        }
        // System.out.println(time() - t_md);
        return new double[] {t_md, loc};
    }
    
    
    public double[] nucleationProbabilityAndLocation(double loc) {
        return new double[] {0, loc};
    }
}
