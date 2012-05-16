package kip.javasim.geometry;

import java.awt.Color;

import scikit.graphics.Drawable;
import scikit.graphics.dim2.Gfx2D;
import scikit.util.Bounds;
import scikit.util.Point;

public class VoronoiGraphics implements Drawable<Gfx2D> {
	private Bounds _bds;
	private QHull _geom;
	private Point[][] _faces;
	
	public VoronoiGraphics(Bounds bds) {
		_bds = bds;
		_geom = new QHull("/sw/bin/qhull");
	}

	public void clear() {
		_faces = null;
	}
	
	public void construct(double[] state, int stride, int N0, int N1) {
		_faces = _geom.constructVoronoi2D(state, stride, N0, N1);
	}
	
	public void draw(Gfx2D g) {
		if (_faces == null)
			return;
		
		g.setColor(Color.RED);
		
		for (Point[] face : _faces) {
			int n = face.length;
			for (int i = 0; i < n; i++) {
				Point v1 = face[(i+0)%n];
				Point v2 = face[(i+1)%n];
				if (v1 != null && v2 != null) {
					g.drawLine(v1.x, v1.y, v2.x, v2.y);
				}
			}
		}
	}
	
	public Bounds getBounds() {
		return _bds;
	}
}
