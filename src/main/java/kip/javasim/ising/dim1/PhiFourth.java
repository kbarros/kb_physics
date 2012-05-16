package kip.javasim.ising.dim1;

import scikit.jobs.params.Parameters;
import static java.lang.Math.*;


public class PhiFourth extends Dynamics1D {
	public double[] field, scratch;
	int N_dx;
	double eps, beta;
    
	public PhiFourth(Parameters params) {
		initialize(params);
	}
	
    public PhiFourth clone() {
		PhiFourth c = (PhiFourth)super.clone();
		c.field = (double[])field.clone();
		return c;
    }
	
	public void initialize(Parameters params) {
		super.initialize(params);
		
		R = params.iget("R");
		dx = (int)(R*params.fget("dx/R"));
		N = (int)(R*params.fget("N/R"));
		N -= N % dx;
		N_dx = N / dx;
		
		field = new double[N_dx];
		scratch = new double[N_dx];
		
		for (int i = 0; i < N_dx; i++)
			field[i] = -1;
	}
	
	public void setParameters(Parameters params) {
		super.setParameters(params);
		eps = params.fget("\u03b5", 1);
		beta = params.fget("beta", 1);
	}
	
	public void randomizeField(double m) {
		setField(m);
	}
	
	public void setField(double m) {
		for (int i = 0; i < N_dx; i++)
			field[i] = m;		
	}
	
	public double magnetization() {
		double sum = 0;
		for (int i = 0; i < N_dx; i++)
			sum += field[i];
		return sum / N_dx;
	}
	
	public double fieldElement(int i) {
		return field[i];
	}

	protected void _step() {
        for (int i = 0; i < N_dx; i++) {
			double phi = field[i];
            double phip = field[(i+1)%N_dx];
            double phim = field[(i-1+N_dx)%N_dx];
            double R_dx = R / (double)dx;
            double R2laplace = R_dx*R_dx*(phim-2*phi+phip);
            double phi3 = phi*phi*phi;
            double eta = random.nextGaussian() * sqrt(2*beta*dt/dx);
			scratch[i] = field[i] + -dt*(-R2laplace + 2*eps*phi + 4*phi3 - h) + eta;
        }
		System.arraycopy(scratch, 0, field, 0, N_dx);
	}
    
    
    // -- NUCLEATION ----------------------------------------------------------------
    //
    
    
    public scikit.dataset.Function saddleProfile() {
        double N_R = (double)N/R;
        return new scikit.dataset.Function(-N_R/2, N_R/2) {
            public double eval(double x) {
                double eps3 = eps*eps*eps;
                double u = 1;
                double h_s = sqrt(-8*eps3/(27*u));
                double phi_s = -sqrt(-eps/(6*u));
                double dh = h_s - h;
                double phi_min = phi_s*(1+sqrt(2*dh/(3*h_s)));
                double a = - sqrt(6*h_s*dh)/(2*phi_s);
                double b = a*a/(3*dh);
                double sch = 1/cosh(x*sqrt(a/2));
                return phi_min + (a/b)*sch*sch;
            }
        };
    }
    
    public double[] nucleationProbabilityAndLocation(double loc) {
        beta = 0; // disable noise
        double t_max = time()+NUCLEATION_WAIT_TIME;
        
        while (!nucleated() && time() < t_max) {
            // Job.plot(0, field);
            step();
        }
        
        beta = 1; // re-enable noise
        if (nucleated())
            return new double[] {1, dx*maxIndex(field)};
        else
            return new double[] {0, Double.NaN};
    }
}
