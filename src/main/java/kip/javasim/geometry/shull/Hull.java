package kip.javasim.geometry.shull;

import java.util.ArrayList;
import java.util.List;

/*
  copyright s-hull.org 2011
  released under the contributors beerware license

  contributors: Phil Atkin, Dr Sinclair.
  translated to Java by Kipton Barros
 */

/// <summary>
/// Vertices belonging to the convex hull need to maintain a point and triad index
/// </summary>
class HullVertex extends Vertex {
	public int pointsIndex;
	public int triadIndex;

	public HullVertex(List<Vertex> points, int pointIndex)
	{
		x = points.get(pointIndex).x;
		y = points.get(pointIndex).y;
		pointsIndex = pointIndex;
		triadIndex = 0;
	}
}

/// <summary>
/// Hull represents a list of vertices in the convex hull, and keeps track of
/// their indices (into the associated points list) and triads
/// </summary>
class Hull extends ArrayList<HullVertex>
{
	private static final long serialVersionUID = -969062917118986409L;

	private int NextIndex(int index)
	{
		if (index == size() - 1)
			return 0;
		else
			return index + 1;
	}

	/// <summary>
	/// Return x component of vector from the hull point at index to next point
	/// </summary>
	public float XDisplacementToNext(int index)
	{
		Vertex et = get(index), en = get(NextIndex(index));
		return en.x - et.x;
	}

	/// <summary>
	/// Return x component of vector from the hull point at index to next point
	/// </summary>
	public float YDisplacementToNext(int index)
	{
		Vertex et = get(index), en = get(NextIndex(index));
		return en.y - et.y;
	}

	/// <summary>
	/// Return whether the hull vertex at index is visible from the supplied coordinates
	/// </summary>
	public boolean EdgeVisibleFrom(int index, float dx, float dy)
	{
		float idx = XDisplacementToNext(index);
		float idy = YDisplacementToNext(index);

		float crossProduct = -dy * idx + dx * idy;
		return crossProduct < 0;
	}

	/// <summary>
	/// Return whether the hull vertex at index is visible from the point
	/// </summary>
	public boolean EdgeVisibleFrom(int index, Vertex point)
	{
		float idx = XDisplacementToNext(index);
		float idy = YDisplacementToNext(index);

		float dx = point.x - get(index).x;
		float dy = point.y - get(index).y;

		float crossProduct = -dy * idx + dx * idy;
		return crossProduct < 0;
	}
}

