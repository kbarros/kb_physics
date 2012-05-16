package kip.javasim.geometry.shull;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class Program {
	public static void main(String[] args) throws NumberFormatException, IOException {
		Random randy = new Random(138);

		List<Vertex> points = new ArrayList<Vertex>();

		if (args.length == 0)
		{
			// Generate random points.
			for (int i = 0; i < 40000; i++)
			{
				float x = randy.nextInt(10000);
				float y = randy.nextInt(10000);
				points.add(new Vertex(x, y));
			}
		}
		else
		{
			// Read a points file as used by the Delaunay Triangulation Tester program DTT
			// (See http://gemma.uni-mb.si/dtt/)
			BufferedReader reader = new BufferedReader(new FileReader(args[0]));
			int count = Integer.parseInt(reader.readLine());
			for (int i = 0; i < count; i++)
			{
				String line = reader.readLine();
				String[] split = line.split("\\s+");
				points.add(new Vertex(Float.parseFloat(split[0]), Float.parseFloat(split[1])));
			}
			reader.close();
		}

		// Write the points in the format suitable for DTT
		PrintWriter writer = new PrintWriter (new FileWriter ("Triangulation-java.pnt"));
		writer.println(points.size());
		for (int i = 0; i < points.size(); i++)
			writer.println(String.format("%f %f", points.get(i).x, points.get(i).y));
		writer.close();
		
		
		// Write out the data set we're actually going to triangulate
		Triangulator angulator = new Triangulator();

		long ms = System.currentTimeMillis();
		List<Triad> triangles = angulator.Triangulation(points, true);
		ms = System.currentTimeMillis() - ms;

		System.out.println("Elapsed time: " + ms + "ms");

		// Write the triangulation results in the format suitable for DTT
		writer = new PrintWriter (new FileWriter ("Triangulation-java.dtt"));
		writer.println(triangles.size());
		for (int i = 0; i < triangles.size(); i++)
		{
			Triad t = triangles.get(i);
			writer.println(String.format("%d: %d %d %d %d %d %d",
					i + 1,
					t.a, t.b, t.c,
					t.ab + 1, t.bc + 1, t.ac + 1));
		}
		writer.close();
	}
}
