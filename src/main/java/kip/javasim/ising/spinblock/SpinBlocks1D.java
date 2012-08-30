package kip.javasim.ising.spinblock;


public class SpinBlocks1D implements Cloneable {
	SpinBlockIndexer indexer;
	int[] indices;
	int[][] blocks;
	int netSum;
	public int L, R;
	
	
    public Object clone() {
        try {
            SpinBlocks1D c = (SpinBlocks1D)super.clone();
            c.indices = (int[])indices.clone(); // not necessary?
            c.blocks = (int[][])blocks.clone();
			for (int i = 0; i < c.blocks.length; i++) {
				c.blocks[i] = (int[])blocks[i].clone();
			}
            return c;
        } catch (Exception e) {
            return null;
        }
    }
	
	
	// L is the system size
	// 2*R+1 is the number of spin in the interaction range.
	// initialize all spins up in direction dir = +-1
	public SpinBlocks1D(int L, int R, int dir) {
		this.L = L;
		this.R = R;
		netSum = L*dir;
		
		indexer = new SpinBlockIndexer(L, R);
		int maxScale = indexer.maxScale();
		indices = indexer.newArray();
		blocks = new int[maxScale+1][];
		
		for (int scale = 0; scale <= maxScale; scale++) {
			int blockLen = 1 << scale;
			blocks[scale] = new int[L/blockLen];
			for (int x = 0; x < L/blockLen; x++) {
				blocks[scale][x] = dir*blockLen; // every block is completely filled in direction dir
			}
		}
	}
	
	
	public int sumAll() {
		return netSum;
	}
	
	public int sumInRange(int x) {
		return sumInRange(x-R, x+R);
	}
	
	
	public int sumInRange(int xlo, int xhi) {
		indexer.fillArray(xlo, xhi, indices);
		int sum = 0;
		for (int i = 0; indices[i] >= 0; i += 2) {
			sum += blocks[indices[i]][indices[i+1]];
		}
		// assert(sum == slowSumInRange(x));
		return sum;
	}
	
	
	public int slowSumInRange(int xlo, int xhi) {
		int sum = 0;
		for (int xp = xlo; xp <= xhi; xp++) {
			int i = (xp + L)%L;
			sum += blocks[0][i];
		}
		return sum;
	}
	
	
	public void flip(int x) {
		int dm = -2*blocks[0][x]; // change in magnetization after flipping spin i
		for (int scale = 0; scale < blocks.length; scale++) {
			int b = x >> scale;
			blocks[scale][b] += dm;
		}
		netSum += dm;
	}
	
	public void set(int x, int s) {
		assert (s == 1 || s == -1);
		if (get(x) != s)
			flip(x);
	}
	
	public int get(int x) {
		return blocks[0][x];
	}
	
	public int getBlock(int scale, int x) {
		return blocks[scale][x];
	}
	
	public int[] getAll() {
		return blocks[0];
	}
}
