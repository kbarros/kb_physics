package kip.javasim.ising.dim2;

//import kip.util.Random;

public class IsingPacked2D extends Ising2D {
//	public int spin[];
	private int packed[];
	
//	public int L1, L2;
//	public int N;
//	public double T;
//	public Random random = new Random();
//	public double time;
//	public boolean openBoundary = false;
	
	int ODD_MASK = 0x55555555;    // = 0b0101...0101
	int EVN_MASK = ODD_MASK << 1; // = 0b1010...1010
	
	
	public IsingPacked2D(int seed, int _L1, int _L2, double _T, boolean _openBoundary) {
		super(seed, _L1, _L2, _T, _openBoundary);
		
		if (L1 % 32 != 0) {
			throw new IllegalArgumentException("System width must be multiple of 32");
		}
		if (_openBoundary) {
			throw new IllegalArgumentException("Periodic boundary conditions only");
		}
		
//		random.setSeed(seed);
//		L1 = _L1;
//		L2 = _L2;
//		T = _T;
//		N = L1*L2;
//		time = 0;
//		spin = new int[N];
		
		packed = new int[N/32];
//		randomize();
		pack();
	}
	
//	public void randomize() {
//		for (int i = 0; i < N/32; i++)
//			packed[i] = random.nextInt();
//		unpack();
//	}
	
	public void pack() {
		for (int i = 0; i < N/32; i++) {
			int x = (i*32) % L1;
			int y = (i*32) / L1;
			
			int a = 0;
			for (int b = 0; b < 32; b++)
				a |= ((spin[y*L1+(x+b)]+1) / 2) << (31-b);
			packed[i] = a;
		}
	}

	public void unpack() {
		for (int i = 0; i < N/32; i++) {
			int x = (i*32) % L1;
			int y = (i*32) / L1;
			
			for (int b = 0; b < 32; b++)
				spin[y*L1+(x+b)] = 2*((packed[i]>>>(31-b))&1)-1;
		}
	}

	
	private int updateEven(int here, int right, int left, int up, int down) {			
		int vert = ((EVN_MASK&up)>>>1) + ((EVN_MASK&down)>>>1);
		int horiz = (ODD_MASK&here) + ((ODD_MASK&here)>>>2) + ((ODD_MASK&left)<<30); 

		int shouldFlip = random.nextInt();
		for (int j = 0; j < 16; j++) {
			if (((shouldFlip>>>(j+15))&1) == 1) {
				int neighborSum = 2*(((vert>>>2*j)&0x3) + ((horiz>>>2*j)&0x3)) - 4;
				int s = 2*((here>>>(2*j+1))&1)-1;
				double dE = 2*s*neighborSum;
				if (dE <= 0 || (T != 0 && random.nextDouble() < Math.exp(-dE/T))) {
					here ^= 1 << (2*j+1);
				}
			}
		}
		return here;
	}
	
	private int updateOdd(int here, int right, int left, int up, int down) {			
			int vert = (ODD_MASK&up) + (ODD_MASK&down);
			int horiz = ((EVN_MASK&here)>>>1) + ((EVN_MASK&here)<<1) + ((EVN_MASK&right)>>>31); 
			
			int shouldFlip = random.nextInt();
			for (int j = 0; j < 16; j++) {
				if (((shouldFlip>>>(j+15))&1) == 1) {
					int neighborSum = 2*(((vert>>>2*j)&0x3) + ((horiz>>>2*j)&0x3)) - 4;
					int s = 2*((here>>>(2*j))&1)-1;
					double dE = 2*s*neighborSum;
					if (dE <= 0 || (T != 0 && random.nextDouble() < Math.exp(-dE/T))) {
						here ^= 1 << (2*j);
					}
				}
			}
			return here;
	}
	
	private void update(int parity) {
		for (int i = 0; i < N/32; i++) {
			int x = i % (L1/32);
			int y = i / (L1/32);
			
			int xp = (x+1+L1/32)%(L1/32);
			int xm = (x-1+L1/32)%(L1/32);
			int yp = (y+1+L2)%(L2);
			int ym = (y-1+L2)%(L2);
			
			int here = packed[i];
			int right = packed[y*(L1/32) + xp];
			int left = packed[y*(L1/32) + xm];
			int up = packed[yp*(L1/32) + x];
			int down = packed[ym*(L1/32) + x];
			
			if (y%2 == parity)
				packed[i] = updateEven(here, right, left, up, down);
			else
				packed[i] = updateOdd(here, right, left, up, down);
		}
	}
	
//	private int neighborSumSlow(int i) {
//		int x = i % L1;
//		int y = i / L1;
//		int acc = 0;
//		
//		if (openBoundary) {
//			if (y < L2-1)
//				acc += spin[i+L1];
//			if (y > 0)
//				acc += spin[i-L1];
//			if (x < L1-1)
//				acc += spin[i+1];
//			if (x > 0)
//				acc += spin[i-1];
//		}
//		else {
//			int yp = (y+1)%L2;
//			int ym = (y-1+L2)%L2;
//			int xp = (x+1)%L1;
//			int xm = (x-1+L1)%L1;
//			acc += spin[yp*L1+x];
//			acc += spin[ym*L1+x];
//			acc += spin[y*L1+xp];
//			acc += spin[y*L1+xm];
//		}
//		
//		return acc;
//	}
//	
//	private void updateSlow(int parity) {
//		for (int i = 0; i < N; i++) {
//			int x = i % L1;
//			int y = i / L1;
//			if ((x + y) % 2 == parity) {
//				double dE = 2*spin[i]*neighborSumSlow(i);
//				if (dE <= 0 || (T != 0 && random.nextDouble() < Math.exp(-dE/T))) {
//					spin[i] = -spin[i];
//				}				
//			}
//		}
//	}
	
	public void step(double mcs) {
		long imcs = Math.round(mcs);
		for (int i = 0; i < imcs; i++) {
			update(0); update(1);
//			updateSlow(0); updateSlow(1); pack();
		}
		time += mcs;
		unpack();
	}
}
