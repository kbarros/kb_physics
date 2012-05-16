package kip.javasim.ising.spinblock;


public class SpinBlocks2D {
	SpinBlockIndexer indexer;
	int[] xIndices, yIndices;
	int[/*yscale*/][/*xscale*/][] blocks;
	int netSum;	
	public int N, L, R;
	
	public SpinBlocks2D(int L, int R) {
		this.L = L;
		this.R = R;
		N = L*L;
		netSum = N;
		
		indexer = new SpinBlockIndexer(L, R);
		int maxScale = indexer.maxScale();
		xIndices = indexer.newArray();
		yIndices = indexer.newArray();
		blocks = new int[maxScale+1][maxScale+1][];
		
		for (int yscale = 0; yscale <= maxScale; yscale++) {
			for (int xscale = 0; xscale <= maxScale; xscale++) {
				int blockSize = 1 << (xscale + yscale);
				blocks[yscale][xscale] = new int[N/blockSize];
				for (int i = 0; i < N/blockSize; i++) {
					blocks[yscale][xscale][i] = blockSize;
				}
			}
		}
	}
	
	
	public int sumAll() {
		return netSum;
	}
	
	public int sumInRange(int x, int y) {
		return sumInRange(x-R, x+R, y-R, y+R);
	}
	
	
	public int sumInRange(int xlo, int xhi, int ylo, int yhi) {
		indexer.fillArray(xlo, xhi, xIndices);
		indexer.fillArray(ylo, yhi, yIndices);
		int sum = 0;
		for (int i = 0; yIndices[i] >= 0; i += 2) {
			int yscale = yIndices[i];
			int yind = yIndices[i+1];
			for (int j = 0; xIndices[j] >= 0; j += 2) {
				int xscale = xIndices[j];
				int xind = xIndices[j+1];
				sum += blocks[yscale][xscale][yind*(L>>xscale) + xind];
			}
		}
		// assert(sum == slowSumInRange(x,y));
		return sum;
	}
	
	
	public int slowSumInRange(int xlo, int xhi, int ylo, int yhi) {
		int sum = 0;
		for (int yp = ylo; yp <= yhi; yp++) {
			for (int xp = xlo; xp <= xhi; xp++) {
				int xi = (xp + L)%L;
				int yi = (yp + L)%L;
				sum += blocks[0][0][yi*L+xi];
			}
		}
		return sum;
	}
	
	
	public void flip(int x, int y) {
		int dm = -2*blocks[0][0][y*L+x]; // change in magnetization after flipping spin i
		for (int yscale = 0; yscale < blocks.length; yscale++) {
			int yind = y >> yscale;
			for (int xscale = 0; xscale < blocks.length; xscale++) {
				int xind = x >> xscale;
				blocks[yscale][xscale][yind*(L>>xscale) + xind] += dm;
			}
		}
		netSum += dm;
	}
	
	public void set(int x, int y, int s) {
		assert (s == 1 || s == -1);
		if (get(x, y) != s)
			flip(x, y);
	}
	
	public int get(int x, int y) {
		return blocks[0][0][y*L+x];
	}
	
	public int[] blocksAtScale(int scale) {
		return blocks[scale][scale];
	}
}

