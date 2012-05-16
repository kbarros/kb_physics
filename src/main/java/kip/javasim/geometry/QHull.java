package kip.javasim.geometry;

import java.io.*;

import scikit.util.Point;


public class QHull {
	private String path;
	private PrintWriter pw;
	private BufferedReader br;
	
	public QHull(String path) {
		this.path = path;
	}
	
	public Point[][] constructVoronoi2D(double[] state, int stride, int N0, int N1) {
		try {
			Process p = Runtime.getRuntime().exec(path + " v o");
			pw = new PrintWriter(p.getOutputStream());
			br = new BufferedReader(new InputStreamReader(p.getInputStream()));
		} catch (IOException e) {
			System.err.println("QHull not available at " + path);
		}

		pw.println(2); // dimension
		pw.println(N1-N0); // number of points
		for (int i = N0; i < N1; i++) {
			pw.println(state[(2*i+0)*stride] + " " + state[(2*i+1)*stride]); // x y coordinates
		}
		pw.close();
		
		try {
			br.readLine(); // dimension
			String[] line = br.readLine().trim().split("\\s+");
			int nvertices = Integer.valueOf(line[0]);
			int nregions = Integer.valueOf(line[1]);
			
			Point vertices[] = new Point[nvertices];
			for (int i = 0; i < nvertices; i++) {
				line = br.readLine().trim().split("\\s+");
				vertices[i] = new Point(Double.valueOf(line[0]), Double.valueOf(line[1]));
			}
			vertices[0] = null; // this is the vertex at "infinity"
			
			Point regions[][] = new Point[nregions][];
			for (int i = 0; i < nregions; i++) {
				line = br.readLine().trim().split("\\s+");
				int npts = line.length-1;
				regions[i] = new Point[npts];
				for (int j = 0; j < npts; j++) {
					regions[i][j] = vertices[Integer.valueOf(line[j+1])];
				}
			}
			return regions;
		} catch (IOException e) {
			System.err.println("Invalid QHull output.");
		}
		return null;
	}
}
