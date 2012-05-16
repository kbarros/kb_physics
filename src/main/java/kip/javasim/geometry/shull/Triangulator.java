package kip.javasim.geometry.shull;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

/*
  copyright s-hull.org 2011
  released under the contributors beerware license

  contributors: Phil Atkin, Dr Sinclair.
 */
class Triangulator
{
	private List<Vertex> points;

	public Triangulator()
	{
	}

	private void Analyse(List<Vertex> suppliedPoints, Hull hull, List<Triad> triads, boolean rejectDuplicatePoints, boolean hullOnly)
	{
		if (suppliedPoints.size() < 3)
			throw new IllegalArgumentException("Number of points supplied must be >= 3");

		this.points = suppliedPoints;
		int nump = points.size();

		float[] distance2ToCentre = new float[nump];
		int[] sortedIndices = new int[nump];

		// Choose first point as the seed
		for (int k = 0; k < nump; k++)
		{
			distance2ToCentre[k] = points.get(0).distance2To(points.get(k));
			sortedIndices[k] = k;
		}

		// Sort by distance to seed point
		QuickSort.quicksort(distance2ToCentre, sortedIndices);

		// Duplicates are more efficiently rejected now we have sorted the vertices
		if (rejectDuplicatePoints)
		{
			// Search backwards so each removal is independent of any other
			for (int k = nump - 2; k >= 0; k--)
			{
				// If the points are identical then their distances will be the same,
				// so they will be adjacent in the sorted list
				if ((points.get(sortedIndices[k]).x == points.get(sortedIndices[k + 1]).x) &&
						(points.get(sortedIndices[k]).y == points.get(sortedIndices[k + 1]).y))
				{
					// Duplicates are expected to be rare, so this is not particularly efficient
					System.arraycopy(sortedIndices, k + 2, sortedIndices, k + 1, nump - k - 2);
					System.arraycopy(distance2ToCentre, k + 2, distance2ToCentre, k + 1, nump - k - 2);
					nump--;
				}
			}
		}

		System.out.println((points.size() - nump) + " duplicate points rejected");

		if (nump < 3)
			throw new IllegalArgumentException("Number of unique points supplied must be >= 3");

		int mid = -1;
		float romin2 = Float.MAX_VALUE, circumCentreX = 0, circumCentreY = 0;

		// Find the point which, with the first two points, creates the triangle with the smallest circumcircle
		Triad tri = new Triad(sortedIndices[0], sortedIndices[1], 2);
		for (int kc = 2; kc < nump; kc++)
		{
			tri.c = sortedIndices[kc];
			if (tri.FindCircumcirclePrecisely(points) && tri.circumcircleR2 < romin2)
			{
				mid = kc;
				// Centre of the circumcentre of the seed triangle
				romin2 = tri.circumcircleR2;
				circumCentreX = tri.circumcircleX;
				circumCentreY = tri.circumcircleY;
			}
			else if (romin2 * 4 < distance2ToCentre[kc])
				break;
		}

		// Change the indices, if necessary, to make the 2th point produce the smallest circumcircle with the 0th and 1th
		if (mid != 2)
		{
			int indexMid = sortedIndices[mid];
			float distance2Mid = distance2ToCentre[mid];

			System.arraycopy(sortedIndices, 2, sortedIndices, 3, mid - 2);
			System.arraycopy(distance2ToCentre, 2, distance2ToCentre, 3, mid - 2);
			sortedIndices[2] = indexMid;
			distance2ToCentre[2] = distance2Mid;
		}

		// These three points are our seed triangle
		tri.c = sortedIndices[2];
		tri.MakeClockwise(points);
		tri.FindCircumcirclePrecisely(points);

		// Add tri as the first triad, and the three points to the convex hull
		triads.add(tri);
		hull.add(new HullVertex(points, tri.a));
		hull.add(new HullVertex(points, tri.b));
		hull.add(new HullVertex(points, tri.c));

		// Sort the remainder according to their distance from its centroid
		// Re-measure the points' distances from the centre of the circumcircle
		Vertex centre = new Vertex(circumCentreX, circumCentreY);
		for (int k = 3; k < nump; k++)
			distance2ToCentre[k] = points.get(sortedIndices[k]).distance2To(centre);

		// Sort the _other_ points in order of distance to circumcentre
		QuickSort.quicksort(distance2ToCentre, sortedIndices, 3, nump - 1);

		// Add new points into hull (removing obscured ones from the chain)
		// and creating triangles....
		int numt = 0;
		for (int k = 3; k < nump; k++)
		{
			int pointsIndex = sortedIndices[k];
			HullVertex ptx = new HullVertex(points, pointsIndex);

			float dx = ptx.x - hull.get(0).x, dy = ptx.y - hull.get(0).y;  // outwards pointing from hull[0] to pt.

			int numh = hull.size();
			List<Integer> pidx = new ArrayList<Integer>(), tridx = new ArrayList<Integer>();
			int hidx;  // new hull point location within hull.....

			if (hull.EdgeVisibleFrom(0, dx, dy))
			{
				// starting with a visible hull facet !!!
				hidx = 0;

				// check to see if segment numh is also visible
				if (hull.EdgeVisibleFrom(numh - 1, dx, dy))
				{
					// visible.
					pidx.add(hull.get(numh - 1).pointsIndex);
					tridx.add(hull.get(numh - 1).triadIndex);

					for (int h = 0; h < numh - 1; h++)
					{
						// if segment h is visible delete h
						pidx.add(hull.get(h).pointsIndex);
						tridx.add(hull.get(h).triadIndex);
						if (hull.EdgeVisibleFrom(h, ptx))
						{
							hull.remove(h);
							h--;
							numh--;
						}
						else
						{
							// quit on invisibility
							hull.add(0, ptx);
							numh++;
							break;
						}
					}
					// look backwards through the hull structure
					for (int h = numh - 2; h > 0; h--)
					{
						// if segment h is visible delete h + 1
						if (hull.EdgeVisibleFrom(h, ptx))
						{
							pidx.add(0, hull.get(h).pointsIndex);
							tridx.add(0, hull.get(h).triadIndex);
							hull.remove(h + 1);  // erase end of chain
						}
						else
							break; // quit on invisibility
					}
				}
				else
				{
					hidx = 1;  // keep pt hull[0]
					tridx.add(hull.get(0).triadIndex);
					pidx.add(hull.get(0).pointsIndex);

					for (int h = 1; h < numh; h++)
					{
						// if segment h is visible delete h  
						pidx.add(hull.get(h).pointsIndex);
						tridx.add(hull.get(h).triadIndex);
						if (hull.EdgeVisibleFrom(h, ptx))
						{                     // visible
							hull.remove(h);
							h--;
							numh--;
						}
						else
						{
							// quit on invisibility
							hull.add(h, ptx);
							break;
						}
					}
				}
			}
			else
			{
				int e1 = -1, e2 = numh;
				for (int h = 1; h < numh; h++)
				{
					if (hull.EdgeVisibleFrom(h, ptx))
					{
						if (e1 < 0)
							e1 = h;  // first visible
					}
					else
					{
						if (e1 > 0)
						{
							// first invisible segment.
							e2 = h;
							break;
						}
					}
				}

				// triangle pidx starts at e1 and ends at e2 (inclusive).	
				if (e2 < numh)
				{
					for (int e = e1; e <= e2; e++)
					{
						pidx.add(hull.get(e).pointsIndex);
						tridx.add(hull.get(e).triadIndex);
					}
				}
				else
				{
					for (int e = e1; e < e2; e++)
					{
						pidx.add(hull.get(e).pointsIndex);
						tridx.add(hull.get(e).triadIndex);   // there are only n-1 triangles from n hull pts.
					}
					pidx.add(hull.get(0).pointsIndex);
				}

				// erase elements e1+1 : e2-1 inclusive.
				if (e1 < e2 - 1)
					hull.subList(e1 + 1, e2).clear();

				// insert ptx at location e1+1.
				hull.add(e1 + 1, ptx);
				hidx = e1 + 1;
			}

			// If we're only computing the hull, we're done with this point
			if (hullOnly)
				continue;

			int a = pointsIndex, T0;

			int npx = pidx.size() - 1;
			numt = triads.size();
			T0 = numt;

			for (int p = 0; p < npx; p++)
			{
				Triad trx = new Triad(a, pidx.get(p), pidx.get(p + 1));
				trx.FindCircumcirclePrecisely(points);

				trx.bc = tridx.get(p);
				if (p > 0)
					trx.ab = numt - 1;
				trx.ac = numt + 1;

				// index back into the triads.
				Triad txx = triads.get(tridx.get(p));
				if ((trx.b == txx.a && trx.c == txx.b) | (trx.b == txx.b && trx.c == txx.a))
					txx.ab = numt;
				else if ((trx.b == txx.a && trx.c == txx.c) | (trx.b == txx.c && trx.c == txx.a))
					txx.ac = numt;
				else if ((trx.b == txx.b && trx.c == txx.c) | (trx.b == txx.c && trx.c == txx.b))
					txx.bc = numt;

				triads.add(trx);
				numt++;
			}
			// Last edge is on the outside
			triads.get(numt - 1).ac = -1;

			hull.get(hidx).triadIndex = numt - 1;
			if (hidx > 0)
				hull.get(hidx - 1).triadIndex = T0;
			else
			{
				numh = hull.size();
				hull.get(numh - 1).triadIndex = T0;
			}
		}
	}

	/// <summary>
	/// Return the convex hull of the supplied points,
	/// Don't check for duplicate points
	/// </summary>
	/// <param name="points">List of 2D vertices</param>
	/// <returns></returns>
	public List<Vertex> ConvexHull(List<Vertex> points)
	{
		return ConvexHull(points, false);
	}

	/// <summary>
	/// Return the convex hull of the supplied points,
	/// Optionally check for duplicate points
	/// </summary>
	/// <param name="points">List of 2D vertices</param>
	/// <param name="rejectDuplicatePoints">Whether to omit duplicated points</param>
	/// <returns></returns>
	public List<Vertex> ConvexHull(List<Vertex> points, boolean rejectDuplicatePoints)
	{
		Hull hull = new Hull();
		List<Triad> triads = new ArrayList<Triad>();

		Analyse(points, hull, triads, rejectDuplicatePoints, true);

		List<Vertex> hullVertices = new ArrayList<Vertex>();
		for (HullVertex hv : hull)
			hullVertices.add(new Vertex(hv.x, hv.y));

		return hullVertices;
	}

	/// <summary>
	/// Return the Delaunay triangulation of the supplied points
	/// Don't check for duplicate points
	/// </summary>
	/// <param name="points">List of 2D vertices</param>
	/// <returns>Triads specifying the triangulation</returns>
	public List<Triad> Triangulation(List<Vertex> points)
	{
		return Triangulation(points, false);
	}

	/// <summary>
	/// Return the Delaunay triangulation of the supplied points
	/// Optionally check for duplicate points
	/// </summary>
	/// <param name="points">List of 2D vertices</param>
	/// <param name="rejectDuplicatePoints">Whether to omit duplicated points</param>
	/// <returns></returns>
	public List<Triad> Triangulation(List<Vertex> points, boolean rejectDuplicatePoints)
	{
		List<Triad> triads = new ArrayList<Triad>();
		Hull hull = new Hull();

		Analyse(points, hull, triads, rejectDuplicatePoints, false);

		// Now, need to flip any pairs of adjacent triangles not satisfying
		// the Delaunay criterion
		int numt = triads.size();
		boolean[] idsA = new boolean[numt];
		boolean[] idsB = new boolean[numt];

		// We maintain a "list" of the triangles we've flipped in order to propogate any
		// consequent changes
		// When the number of changes is large, this is best maintained as a vector of bools
		// When the number becomes small, it's best maintained as a set
		// We switch between these regimes as the number flipped decreases
		int flipped = FlipTriangles(triads, idsA);

		int iterations = 1;
		while (flipped > (int)(fraction * (float)numt))
		{
			if ((iterations & 1) == 1)
				flipped = FlipTriangles(triads, idsA, idsB);
			else
				flipped = FlipTriangles(triads, idsB, idsA);

			iterations++;
		}

		SortedSet<Integer> idSetA = new TreeSet<Integer>(), idSetB = new TreeSet<Integer>();
		flipped = FlipTriangles(triads,
				((iterations & 1) == 1) ? idsA : idsB, idSetA);

		iterations = 1;
		while (flipped > 0)
		{
			if ((iterations & 1) == 1)
				flipped = FlipTriangles(triads, idSetA, idSetB);
			else
				flipped = FlipTriangles(triads, idSetB, idSetA);

			iterations++;
		}

		return triads;
	}

	public float fraction = 0.3f;

	/// <summary>
	/// Test the triad against its 3 neighbours and flip it with any neighbour whose opposite point
	/// is inside the circumcircle of the triad
	/// </summary>
	/// <param name="triads">The triads</param>
	/// <param name="triadIndexToTest">The index of the triad to test</param>
	/// <param name="triadIndexFlipped">Index of adjacent triangle it was flipped with (if any)</param>
	/// <returns>true iff the triad was flipped with any of its neighbours</returns>
	boolean FlipTriangle(List<Triad> triads, int triadIndexToTest, IntWrapper triadIndexFlipped)
	{
		IntWrapper oppositeVertexW = new IntWrapper();
		int edge1, edge2;
		IntWrapper edge3W = new IntWrapper();
		IntWrapper edge4W = new IntWrapper();

		triadIndexFlipped.value = 0;

		Triad tri = triads.get(triadIndexToTest);
		// test all 3 neighbours of tri 

		if (tri.bc >= 0)
		{
			triadIndexFlipped.value = tri.bc;
			Triad t2 = triads.get(triadIndexFlipped.value);
			// find relative orientation (shared limb).
			t2.FindAdjacency(tri.b, triadIndexToTest, oppositeVertexW, edge3W, edge4W);
			if (tri.InsideCircumcircle(points.get(oppositeVertexW.value)))
			{  // not valid in the Delaunay sense.
				edge1 = tri.ab;
				edge2 = tri.ac;
				if (edge1 != edge3W.value && edge2 != edge4W.value)
				{
					int tria = tri.a, trib = tri.b, tric = tri.c;
					tri.Initialize(tria, trib, oppositeVertexW.value, edge1, edge3W.value, triadIndexFlipped.value, points);
					t2.Initialize(tria, tric, oppositeVertexW.value, edge2, edge4W.value, triadIndexToTest, points);

					// change knock on triangle labels.
					if (edge3W.value >= 0)
						triads.get(edge3W.value).ChangeAdjacentIndex(triadIndexFlipped.value, triadIndexToTest);
					if (edge2 >= 0)
						triads.get(edge2).ChangeAdjacentIndex(triadIndexToTest, triadIndexFlipped.value);
					return true;
				}
			}
		}


		if (tri.ab >= 0)
		{
			triadIndexFlipped.value = tri.ab;
			Triad t2 = triads.get(triadIndexFlipped.value);
			// find relative orientation (shared limb).
			t2.FindAdjacency(tri.a, triadIndexToTest, oppositeVertexW, edge3W, edge4W);
			if (tri.InsideCircumcircle(points.get(oppositeVertexW.value)))
			{  // not valid in the Delaunay sense.
				edge1 = tri.ac;
				edge2 = tri.bc;
				if (edge1 != edge3W.value && edge2 != edge4W.value)
				{
					int tria = tri.a, trib = tri.b, tric = tri.c;
					tri.Initialize(tric, tria, oppositeVertexW.value, edge1, edge3W.value, triadIndexFlipped.value, points);
					t2.Initialize(tric, trib, oppositeVertexW.value, edge2, edge4W.value, triadIndexToTest, points);

					// change knock on triangle labels.
					if (edge3W.value >= 0)
						triads.get(edge3W.value).ChangeAdjacentIndex(triadIndexFlipped.value, triadIndexToTest);
					if (edge2 >= 0)
						triads.get(edge2).ChangeAdjacentIndex(triadIndexToTest, triadIndexFlipped.value);
					return true;
				}
			}
		}

		if (tri.ac >= 0)
		{
			triadIndexFlipped.value = tri.ac;
			Triad t2 = triads.get(triadIndexFlipped.value);
			// find relative orientation (shared limb).
			t2.FindAdjacency(tri.a, triadIndexToTest, oppositeVertexW, edge3W, edge4W);
			if (tri.InsideCircumcircle(points.get(oppositeVertexW.value)))
			{  // not valid in the Delaunay sense.
				edge1 = tri.ab;   // .ac shared limb
				edge2 = tri.bc;
				if (edge1 != edge3W.value && edge2 != edge4W.value)
				{
					int tria = tri.a, trib = tri.b, tric = tri.c;
					tri.Initialize(trib, tria, oppositeVertexW.value, edge1, edge3W.value, triadIndexFlipped.value, points);
					t2.Initialize(trib, tric, oppositeVertexW.value, edge2, edge4W.value, triadIndexToTest, points);

					// change knock on triangle labels.
					if (edge3W.value >= 0)
						triads.get(edge3W.value).ChangeAdjacentIndex(triadIndexFlipped.value, triadIndexToTest);
					if (edge2 >= 0)
						triads.get(edge2).ChangeAdjacentIndex(triadIndexToTest, triadIndexFlipped.value);
					return true;
				}
			}
		}

		return false;
	}

	/// <summary>
	/// Flip triangles that do not satisfy the Delaunay condition
	/// </summary>
	private int FlipTriangles(List<Triad> triads, boolean[] idsFlipped)
	{
		int numt = triads.size();
		for (int i = 0; i < numt; i++)
			idsFlipped[i] = false;

		int flipped = 0;
		for (int t = 0; t < numt; t++)
		{
			IntWrapper t2W = new IntWrapper();
			if (FlipTriangle(triads, t, t2W))
			{
				flipped += 2;
				idsFlipped[t] = true;
				idsFlipped[t2W.value] = true;

			}
		}

		return flipped;
	}

	private int FlipTriangles(List<Triad> triads, boolean[] idsToTest, boolean[] idsFlipped)
	{
		int numt = triads.size();
		for (int i = 0; i < numt; i++)
			idsFlipped[i] = false;

		int flipped = 0;
		for (int t = 0; t < numt; t++)
		{
			if (idsToTest[t])
			{
				IntWrapper t2W = new IntWrapper();
				if (FlipTriangle(triads, t, t2W))
				{
					flipped += 2;
					idsFlipped[t] = true;
					idsFlipped[t2W.value] = true;
				}
			}
		}

		return flipped;
	}

	private int FlipTriangles(List<Triad> triads, boolean[] idsToTest, SortedSet<Integer> idsFlipped)
	{
		int numt = triads.size();
		idsFlipped.clear();

		int flipped = 0;
		for (int t = 0; t < numt; t++)
		{
			if (idsToTest[t])
			{
				IntWrapper t2W = new IntWrapper();
				if (FlipTriangle(triads, t, t2W))
				{
					flipped += 2;
					idsFlipped.add(t);
					idsFlipped.add(t2W.value);
				}
			}
		}

		return flipped;
	}

	private int FlipTriangles(List<Triad> triads, SortedSet<Integer> idsToTest, SortedSet<Integer> idsFlipped)
	{
		int flipped = 0;
		idsFlipped.clear();

		for (int t : idsToTest)
		{
			IntWrapper t2W = new IntWrapper();
			if (FlipTriangle(triads, t, t2W))
			{
				flipped += 2;
				idsFlipped.add(t);
				idsFlipped.add(t2W.value);
			}
		}

		return flipped;
	}

	//
	// methods to verify that triad adjacency and indices are set correctly
	//

	void VerifyHullContains(Hull hull, int idA, int idB)
	{
		if (
				((hull.get(0).pointsIndex == idA) && (hull.get(hull.size() - 1).pointsIndex == idB)) ||
				((hull.get(0).pointsIndex == idB) && (hull.get(hull.size() - 1).pointsIndex == idA)))
			return;

		for (int h = 0; h < hull.size() - 1; h++)
		{
			if (hull.get(h).pointsIndex == idA)
			{
				assert(hull.get(h + 1).pointsIndex == idB);
				return;
			}
			else if (hull.get(h).pointsIndex == idB)
			{
				assert(hull.get(h + 1).pointsIndex == idA);
				return;
			}
		}

	}

	void VerifyTriadContains(Triad tri, int nbourTriad, int idA, int idB)
	{
		if (tri.ab == nbourTriad)
		{
			assert(
					((tri.a == idA) && (tri.b == idB)) ||
					((tri.b == idA) && (tri.a == idB)));
		}
		else if (tri.ac == nbourTriad)
		{
			assert(
					((tri.a == idA) && (tri.c == idB)) ||
					((tri.c == idA) && (tri.a == idB)));
		}
		else if (tri.bc == nbourTriad)
		{
			assert(
					((tri.c == idA) && (tri.b == idB)) ||
					((tri.b == idA) && (tri.c == idB)));
		}
		else
			assert(false);
	}

	void VerifyTriads(List<Triad> triads, Hull hull)
	{
		for (int t = 0; t < triads.size(); t++)
		{
			if (t == 17840)
				t = t + 0;

			Triad tri = triads.get(t);
			if (tri.ac == -1)
				VerifyHullContains(hull, tri.a, tri.c);
			else
				VerifyTriadContains(triads.get(tri.ac), t, tri.a, tri.c);

			if (tri.ab == -1)
				VerifyHullContains(hull, tri.a, tri.b);
			else
				VerifyTriadContains(triads.get(tri.ab), t, tri.a, tri.b);

			if (tri.bc == -1)
				VerifyHullContains(hull, tri.b, tri.c);
			else
				VerifyTriadContains(triads.get(tri.bc), t, tri.b, tri.c);

		}
	}

	void WriteTriangles(List<Triad> triangles, String name) throws IOException
	{
		PrintWriter writer = new PrintWriter(new FileWriter(name));
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

	static class QuickSort {
		public static void quicksort(float[] main, int[] index) {
			quicksort(main, index, 0, index.length - 1);
		}

		// quicksort a[left] to a[right]
		public static void quicksort(float[] a, int[] index, int left, int right) {
			if (right <= left) return;
			int i = partition(a, index, left, right);
			quicksort(a, index, left, i-1);
			quicksort(a, index, i+1, right);
			
//			for (int j = left; j < right; j++) {
//				assert(a[j] <= a[j+1]);
//			}
		}

		// partition a[left] to a[right], assumes left < right
		private static int partition(float[] a, int[] index, 
				int left, int right) {
			int i = left - 1;
			int j = right;
			while (true) {
				while (less(a[++i], a[right]))      // find item on left to swap
					;                               // a[right] acts as sentinel
				while (less(a[right], a[--j]))      // find item on right to swap
					if (j == left) break;           // don't go out-of-bounds
				if (i >= j) break;                  // check if pointers cross
				exch(a, index, i, j);               // swap two elements into place
			}
			exch(a, index, i, right);               // swap with partition element
			return i;
		}

		// is x < y ?
		private static boolean less(float x, float y) {
			return (x < y);
		}

		// exchange a[i] and a[j]
		private static void exch(float[] a, int[] index, int i, int j) {
			float swap = a[i];
			a[i] = a[j];
			a[j] = swap;
			int b = index[i];
			index[i] = index[j];
			index[j] = b;
		}
	}
}
