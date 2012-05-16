package kip.javasim.electro.apps;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

import scikit.dataset.Accumulator;
import scikit.util.Commands;
import static java.lang.Math.*; 

class Charges {
	// position of charge
	double x[];
	double y[];
	double z[];
	
	// amount of charge
	double q[];

	Charges(int n) {
		x = new double[n];
		y = new double[n];
		z = new double[n];
		q = new double[n];
	}
}


class SurfaceCharges extends Charges {
	// surface normal unit vector
	double nx[];
	double ny[];
	double nz[];
	
	// surface patch area
	double a[];
	
	// relative permittivity inside and outside of surface
	double e_in[];
	double e_out[];

	SurfaceCharges(int n) {
		super(n);
		nx = new double[n];
		ny = new double[n];
		nz = new double[n];
		a = new double[n];
		e_in = new double[n];
		e_out = new double[n];
	}
	
	void copyArrays(double[] ret, double[] a, double[] b) {
		for (int i = 0; i < a.length; i++)
			ret[i] = a[i];
		for (int i = 0; i < b.length; i++)
			ret[a.length+i] = b[i];		
	}
	
	SurfaceCharges add(SurfaceCharges c) {
		SurfaceCharges ret = new SurfaceCharges(x.length + c.x.length);
		
		copyArrays(ret.x, x, c.x);
		copyArrays(ret.y, y, c.y);
		copyArrays(ret.z, z, c.z);
		copyArrays(ret.q, q, c.q);
		
		copyArrays(ret.nx, nx, c.nx);
		copyArrays(ret.ny, ny, c.ny);
		copyArrays(ret.nz, nz, c.nz);
		copyArrays(ret.a, a, c.a);
		copyArrays(ret.e_in, e_in, c.e_in);
		copyArrays(ret.e_out, e_out, c.e_out);
		
		return ret;
	}
	
	double totalCharge() {
		double net_q = 0;
		for (int i = 0; i < x.length; i++) {
			net_q += q[i];
		}
		return net_q;
	}
}


public class FHMApp {
	
	double[] electricFieldAtPoint(double x, double y, double z, Charges... cs) {
		double[] E = new double[]{0, 0, 0};
		
		for (Charges c : cs) {
			if (c != null) {
				for (int i = 0; i < c.x.length; i++) {
					double dx = x - c.x[i];
					double dy = y - c.y[i];
					double dz = z - c.z[i];

					double r2 = dx*dx + dy*dy + dz*dz;
					double r3 = r2 * Math.sqrt(r2);
					if (r2 == 0) continue; // avoid singularity 

					E[0] += (c.q[i] / (4. * Math.PI)) * dx / r3;
					E[1] += (c.q[i] / (4. * Math.PI)) * dy / r3;
					E[2] += (c.q[i] / (4. * Math.PI)) * dz / r3;
				}
			}
		}
		
		return E;
	}
	
	SurfaceCharges createSurfaceChargesForSlab(double dx, int lp, double e_in) {
		int cnt = 2*lp*lp;
		double e_out = 1.;
		
		SurfaceCharges h = new SurfaceCharges(cnt);
		
		int i = 0;
		
		// create plates with corners at (+-dx/2,-1/2,-1/2) and (+-dx/2, 1/2, 1/2)
		for (int sign = -1; sign <= 1; sign += 2) {
			for (int iy = 0; iy < lp; iy++) {
				for (int iz = 0; iz < lp; iz++) {
					h.x[i] = sign*dx/2.;
					h.y[i] = iy / (lp - 1.0) - 1./2;
					h.z[i] = iz / (lp - 1.0) - 1./2;
					h.q[i] = 0; // surface charge currently unknown

					h.nx[i] = sign;
					h.ny[i] = 0;
					h.nz[i] = 0;
					h.a[i] = (1.*1.) / (lp*lp);
					h.e_in[i] = e_in;
					h.e_out[i] = e_out;

					i++;
				}
			}
		}

		return h;
	}
	
	ArrayList<Double> loadFloats(File f) {
		try {
			ArrayList<Double> ret = new ArrayList<Double>();
			Scanner s = new Scanner(f);
			while (s.hasNextDouble()) {
				ret.add(s.nextDouble());
			}
			return ret;
		} catch (IOException e) {
			System.out.println(e.toString());
			return null;
		}
	}
	
	Charges createPointCharge(double q, double x0, double y0, double z0) {
		Charges ret = new Charges(1);
		ret.q[0] = q;
		ret.x[0] = x0;
		ret.y[0] = y0;
		ret.z[0] = z0;
		return ret;
	}
	
	SurfaceCharges createSurfaceChargesForSphere(int cnt, double e_in, double x0, double y0, double z0) {
		double e_out = 1.;		
		SurfaceCharges h = new SurfaceCharges(cnt);
		
		ArrayList<Double> coords = loadFloats(new File("/Users/kbarros/Desktop/spheres/packing"+cnt+".txt"));
		if (coords.size() != 3*cnt) {
			System.err.println("File sizes don't match");
		}
		
		for (int i = 0; i < cnt; i++) {
			double nx = coords.get(3*i+0);
			double ny = coords.get(3*i+1);
			double nz = coords.get(3*i+2); 
			double r = Math.sqrt(nx*nx+ny*ny+nz*nz);
			nx /= r;
			ny /= r;
			nz /= r;

			h.x[i] = nx + x0;
			h.y[i] = ny + y0;
			h.z[i] = nz + z0;
			h.q[i] = 0; // surface charge currently unknown
			
			h.nx[i] = nx;
			h.ny[i] = ny;
			h.nz[i] = nz;
			
			h.a[i] = 4*Math.PI / cnt;
			h.e_in[i] = e_in;
			h.e_out[i] = e_out;
		}
		
		return h;
	}
	
	SurfaceCharges createSurfaceChargesForCylinder(int cntRequested, double e_in, double x0, double y0, double z0, double len, double rad) {
		double e_out = 1.;
		
		double circ = 2*PI*rad;
		int cnt_len  = (int) sqrt((len/circ)*cntRequested);
		int cnt_circ = (int) sqrt((circ/len)*cntRequested);
		
		SurfaceCharges h = new SurfaceCharges(cnt_len*cnt_circ);
		
		double patch_area = circ*len / (cnt_len*cnt_circ);

		for (int j1 = 0; j1 < cnt_len; j1++) {
			for (int j2 = 0; j2 < cnt_circ; j2++) {
				int i = j1*cnt_circ + j2;
				
		        double theta = ((float)j2/cnt_circ) * 2*PI;
		        
		        double x = rad*cos(theta);
		        double y = rad*sin(theta);
		        double z = ((float)j1/cnt_len) * len;
		        
		        h.x[i] = x + x0;
		        h.y[i] = y + y0;
		        h.z[i] = z + z0;
		        h.q[i] = 0; // surface charge currently unknown
		        // System.out.println(h.x[i] + " " + h.y[i] + " " + h.z[i]);
		        
				// normal vector of cylinder points inwards, indicating that "e_in" actually refers
				// to the volume dielectric, while "e_out" refers to the dielectric in the tube.
		        h.nx[i] = -cos(theta);
		        h.ny[i] = -sin(theta);
		        h.nz[i] = 0;
		        
		        h.a[i] = patch_area;
		        h.e_in[i] = e_in;
		        h.e_out[i] = e_out;
			}
		}
		
		return h;
	}

	void surfaceChargeGradient(double[] grad, SurfaceCharges h, Charges g, double[] E_applied) {
		for (int i = 0; i < h.x.length; i++) {
			double delta_e = h.e_out[i] - h.e_in[i];
			double mean_e = (h.e_out[i] + h.e_in[i]) / 2;
			
			double[] E = electricFieldAtPoint(h.x[i], h.y[i], h.z[i], h, g);
			E[0] += E_applied[0];
			E[1] += E_applied[1];
			E[2] += E_applied[2];
			
			double E_dot_n = E[0]*h.nx[i] + E[1]*h.ny[i] + E[2]*h.nz[i];
			
			grad[i] = h.a[i]*delta_e*E_dot_n + mean_e*h.q[i];
		}
	}
	
	
	Accumulator getChargesByAngle(Charges c) {
		Accumulator ret = new Accumulator();
		
		int numBins = 40;
		double dtheta = Math.PI/numBins;
		for (int i = 0; i < numBins; i++) {
			double theta1 = i*dtheta;
			double theta2 = (i+1)*dtheta;
			double areaOfStrip = 2*Math.PI*(Math.cos(theta1)-Math.cos(theta2));
			
			double charge = 0;
			for (int j = 0; j < c.x.length; j++) {
				double theta = Math.acos(-c.x[j]);
				if (theta1 <= theta && theta < theta2)
					charge += c.q[j];
			}
			
			ret.accum((theta1+theta2)/2., charge/areaOfStrip);
		}
		
		return ret;
	}
	
	
	Accumulator getChargesByZ(SurfaceCharges c, double zmin, double zmax) {
		Accumulator ret = new Accumulator();
		for (int j = 0; j < c.x.length; j++) {
			ret.accum(c.z[j], c.q[j]/c.a[j]);
		}
		return ret;
	}

	
	Accumulator getEField(Charges c, double xlo, double xhi) {
		Accumulator ret = new Accumulator();
		
		int numBins = 40;
		for (int i = 0; i < numBins; i++) {
			double alpha = ((double)i / (numBins-1.0));
			double x = xlo + (xhi-xlo)*alpha;
			ret.accum(x, electricFieldAtPoint(x, 0, 0, c)[0]);
		}
		
		return ret;
	}
	
	public void go() {
		SurfaceCharges h = null;
		Charges g = null;
		double e_in = 2;
		double dt = 0;
		double[] E_applied = null;
		
		int chargeType = 2;
		if (chargeType == 0) {
			dt = 1;
			E_applied = new double[] {1, 0, 0};
			double dx = 0.1;
			int lp = 40;
			h = createSurfaceChargesForSlab(dx, lp, e_in);
		}
		else if (chargeType == 1) {
			dt = 1;
			E_applied = new double[] {1, 0, 0};
//			int cnt = 20; // dodecahedron
//			int cnt = 124; // putatively optimal
//			int cnt = 1002; // icosahedral
//			int cnt = 1592; // putatively optimal
			int cnt = 2012; // icosahedral
			h = createSurfaceChargesForSphere(cnt, e_in, 0, 0, 0);
//			h = h.add(createSurfaceChargesForSphere(cnt, e_in, 2, 2, 0));
		}
		else if (chargeType == 2) {
			dt = 0.5;
			g = createPointCharge(1, 0, 0, 0);
			E_applied = new double[] {0, 0, 0};
			int cnt = 2000;
			double len = 20;
			double x0=0, y0=0, z0=-len/2;
			double rad = 1;
			
			// normal vector of cylinder points inwards, indicating that "e_in" actually refers
			// to the volume dielectric, while "e_out" refers to the dielectric in the tube.
			h = createSurfaceChargesForCylinder(cnt, e_in, x0, y0, z0, len, rad);
		}
		
		System.out.println("Total atoms " + h.x.length);
		
		int nh = h.x.length;
		double[] grad = new double[nh];
		
		for (int cnt = 0; cnt < 10; cnt++) {
			surfaceChargeGradient(grad, h, g, E_applied);
			
//			Commands.plot(getChargesByAngle(h));
			Commands.plot(getChargesByZ(h, -10, 10));
//			Commands.plot(getEField(h, -1.5, 1.5));

			System.out.println(electricFieldAtPoint(0, 0.0, 0, h)[0]);
//			System.out.println(electricFieldAtPoint(-0.5, 0.0, 0, h)[0]);
//			System.out.println(electricFieldAtPoint(-0.9, 0.0, 0, h)[0]);
			
			double max_abs_dq = 0;
			double max_abs_q  = 0;
			double net_dq = 0;
			
			for (int i = 0; i < nh; i++) {
				double dq = - dt * grad[i]; 
				h.q[i] += dq;
				net_dq += dq; 
				
				max_abs_dq = max(max_abs_dq, abs(dq));
				max_abs_q  = max(max_abs_q,  abs(h.q[i]));
			}
			
			
			System.out.println("convergence: " + (max_abs_dq / max_abs_q));
			System.out.println("net_dq: " + net_dq);
			System.out.println("induced charge: " + h.totalCharge());
			System.out.println();
		}
	}
	
	
	public static void main(String[] args) {
		new FHMApp().go();
	}
}
