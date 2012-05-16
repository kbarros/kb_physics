package kip.javasim.md;

import static java.lang.Math.*;

import java.util.ArrayList;

import kip.javasim.Vec3J;
import scikit.numerics.ode.*;


public class MolecularDynamics2D<Pt extends Particle> {
	public ParticleContext pc;
	public Pt[] particles;
	public int N;
	private double time;
	
	protected boolean canonicalEnsemble = false; // microcanonical ensemble by default
	protected double T, gamma; // temperature, and thermodynamic coupling of system to heat bath
	
	
	// complete state of configuration.  positions, velocities, packed as
	// [x_1, vx_1, y_1, vy_1, ..., time]
	protected double phase[];
	// differential equation solver
	protected ODESolver solver;
	// grid for finding neighbors quickly
	JPointGrid2D<Pt> grid;
	
	
	public MolecularDynamics2D(double dt, ParticleContext pc, Pt[] particles) {
		this.pc = pc;
		this.particles = particles;
		N = particles.length;
		time = 0;
		
		// initialize phase space array and ODE solver
		phase = new double[4*N];
		ODE ode = new ODE() {
			public void getRate(double[] state, double[] rate) {
				MolecularDynamics2D.this.getRate(state, rate);
			}
			public double[] getState() {
				return phase;
			}
		};
		grid = new JPointGrid2D<Pt>(pc.L, (int)sqrt(N), pc.periodic(), particles);
		writeStateArray(phase);
		solver = new scikit.numerics.ode.Verlet(ode, 4*N);
		solver.initialize(dt);
	}
	
	public void setStepSize(double dt) {
		solver.setStepSize(dt);
	}
	
	public double getStepSize() {
		return solver.getStepSize();
	}
	
	private void writeStateArray(double[] state) {
		for (int i = 0; i < N; i++) {
			state[4*i+0] = particles[i].x; 
			state[4*i+1] = particles[i].vx;
			state[4*i+2] = particles[i].y;
			state[4*i+3] = particles[i].vy; 
		}
	}
	
	private void readStateArray(double[] state) {
		for (int i = 0; i < N; i++) {
			Pt p = particles[i];
			p.x  = state[4*i+0];
			p.vx = state[4*i+1];
			p.y  = state[4*i+2];
			p.vy = state[4*i+3];
			pc.wrap(p);	
		}
	}
	
	public void step() {
		writeStateArray(phase);
		solver.step();
		readStateArray(phase);
		if (canonicalEnsemble)
			brownianNoise();
		sanityCheck();
		time += solver.getStepSize();
	}
	
	
	public void setTemperature(double T, double gamma) {
		if (T >= 0) {
			canonicalEnsemble = true;
			this.T = T;
			this.gamma = gamma;
		}
		else {
			disableTemperature();
		}
	}
	
	public void disableTemperature() {
		canonicalEnsemble = false;
	}
	
	public double time() {
		return time;
	}
	
	public double potentialEnergy() {
		double V = 0;
		
		grid.initialize();
		for (Pt p1 : particles) {
			for (Pt p2 : grid.pointOffsetsWithinRange(p1, p1.tag.interactionRange)) {
				if (p1 != p2)
					V += p1.potential(p2)/2.; // divisor of 2 corrects for double counting
			}
			// accumulate accelerations due to external forces
			V += p1.potential();
		}
		return V;
	}
	
	
	public double kineticEnergy() {
		double K = 0;
		for (Pt p : particles) {
			double M = p.tag.mass;
			K += 0.5*M*(p.vx*p.vx+p.vy*p.vy);
		}
		return K;
	}
	
	public double reducedKineticEnergy() {
		return (kineticEnergy() - N*T) / N*T;
	}

	
	private void getRate(double[] state, double[] rate) {
		readStateArray(state);
		for (int i = 0; i < N; i++) {
			// set dx/dt = v
			rate[4*i+0] = particles[i].vx;
			rate[4*i+2] = particles[i].vy;
		}
		calculateForces(rate);
	}
	
	private void sanityCheck() {
		for (int i = 0; i < N; i++) {
			Pt p = particles[i];
			double dt = solver.getStepSize();
			if (dt*max(abs(p.vx), abs(p.vy)) > pc.L/10) {
				throw new IllegalStateException("Simulation has destablized");
			}
		}
	}
	
	private void calculateForces(double[] rate) {
		Vec3J f = new Vec3J(0,0,0);
		
		grid.initialize();
		for (int i = 0; i < N; i++) {
			Pt p1 = particles[i];
			double M = p1.tag.mass;
			
			// initialize accelerations to zero
			rate[4*i+1] = 0;
			rate[4*i+3] = 0;
			
			// accumulate accelerations due to pairwise interactions
			ArrayList<Pt> pts = grid.pointOffsetsWithinRange(p1, p1.tag.interactionRange);
			for (int j = 0; j < pts.size(); j++) {
				Pt p2 = pts.get(j);
				if (p1 != p2) {
					p1.force(p2, f);
					rate[4*i+1] += f.x/M;
					rate[4*i+3] += f.y/M;
				}
			}
			
			// accumulate accelerations due to external forces
			p1.force(f);
			rate[4*i+1] += f.x/M;
			rate[4*i+3] += f.y/M;
		}
	}
	
	private void brownianNoise() {
		// accumulate accelerations due to stochastic noise
		for (Pt p : particles) {
			double M = p.tag.mass;
			// dp/dt = - gamma 2p/2m + sqrt(2 gamma T) eta
			// dv/dt = - (gamma v + sqrt(2 gamma T) eta) / m
			double dt = solver.getStepSize();
			p.vx += (-dt*gamma*p.vx + sqrt(2*dt*gamma*T)*pc.rand.nextGaussian())/M;
			p.vy += (-dt*gamma*p.vy + sqrt(2*dt*gamma*T)*pc.rand.nextGaussian())/M;
		}
	}
}
