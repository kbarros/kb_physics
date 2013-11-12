package kip.javasim.ising.dim2.apps;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Scanner;

import kip.javasim.ising.PercolationSite2d;
import scikit.util.Parsing;

public class QuenchedAnalysisApp {
	
	public static int countDegenerates(int w, int h, int[] spin) {
		int cnt = 0;
		for (int y = 0; y < h; y++) {
			for (int x = 0; x < w; x++) {
				int xp = (x+1)%w;
				int yp = (y+1)%h;
				
				if (spin[y*w+x] == spin[yp*w+xp] &&
					spin[yp*w+x] == spin[y*w+xp] &&
					spin[y*w+x] != spin[y*w+xp])
					cnt++;
			}
		}
		return cnt;
	}
	
	public static void main(String[] args) {
		BufferedReader stdin = new BufferedReader(new InputStreamReader(System.in));
		
		int numDegenerates = 0;
		int numGrids = 0;
		int numCross = 0, numPoint = 0, numVert = 0, numHoriz = 0, numDiag = 0;
		
		try {
			String line;
			while ((line = stdin.readLine()) != null) {
				if (line.startsWith("#") || line.trim().equals(""))
					continue;
				
				Scanner sc = new Scanner(line);
				int w = sc.nextInt();
				int h = sc.nextInt();
				assert(!sc.hasNextInt());
				sc.close();
				
				int[] spin = new int[w*h];
				PercolationSite2d nz = new PercolationSite2d(w, h, false); // periodic boundary conditions

				String[] vs = new String[w];
				for (int y = 0; y < h; y++) {
					line = stdin.readLine();
					int splitCnt = Parsing.stringSplit(vs, line);
					assert(splitCnt == w);
					
					for (int x = 0; x < w; x++) {
						if (vs[x].equals("1")) {
							spin[y*w+x] = 1;
						}
						else if (vs[x].equals("0")) {
							spin[y*w+x] = 0;
						}
						else {
							System.out.println("Invalid spin entry: " + vs[x]);
						}
					}
					
				}
				
				numGrids++;
				numDegenerates += countDegenerates(w, h, spin);
				
				nz.occupyAndBondSites(spin, 1);
				nz.findHomologies();
				if (nz.crossHomology()) {
					System.out.println("cross");
					numCross++;
				}
				else if (nz.pointHomology()) {
					System.out.println("point");
					numPoint++;
				}
				else if (nz.verticalHomology()) {
					System.out.println("vert");
					numVert++;
				}
				else if (nz.horizontalHomology()) {
					System.out.println("horiz");
					numHoriz++;
				}
				else {
					System.out.println("diag");
					numDiag++;
				}
			}
			
			System.err.println(
						"# Analysis of "+numGrids+" grids completed.\n"+
						"# Found "+numDegenerates+" degenerate points where four domains meet.\n"+
						"# Summary: Cross="+numCross+" Point="+numPoint+" Vert="+numVert+" Horiz="+numHoriz+" Diag="+numDiag
			);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
