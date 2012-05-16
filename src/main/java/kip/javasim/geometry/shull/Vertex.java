package kip.javasim.geometry.shull;


/*
  copyright s-hull.org 2011
  released under the contributors beerware license

  contributors: Phil Atkin, Dr Sinclair.
  translated to Java by Kipton Barros
*/

public class Vertex
{
	public float x, y;

	protected Vertex() { }

	public Vertex(float x, float y) 
	{
		this.x = x; this.y = y;
	}

	public float distance2To(Vertex other)
	{
		float dx = x - other.x;
		float dy = y - other.y;
		return dx * dx + dy * dy;
	}

	public float distanceTo(Vertex other)
	{
		return (float)Math.sqrt(distance2To(other));
	}

	public String toString()
	{
		return String.format("(%f,%f)", x, y);
	}
}
