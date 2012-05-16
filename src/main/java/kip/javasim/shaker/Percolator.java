package kip.javasim.shaker;

import java.util.ArrayList;

public class Percolator<T> {
	ArrayList<ArrayList<T>> groups;
	
	public Percolator() {
		groups = new ArrayList<ArrayList<T>>();
	}
	
	public void add(T obj) {
		ArrayList<T> g = new ArrayList<T>();
		g.add(obj);
		groups.add(g);
	}
	
	public void bond(T obj1, T obj2) {
		ArrayList<T> g1=null, g2=null;
		for (ArrayList<T> g : groups) {
			if (g.contains(obj1))
				g1 = g;
			if (g.contains(obj2))
				g2 = g;
		}
		if (g1 == null)
			throw new IllegalArgumentException("Object " + obj1 + " not contained in Percolator.");
		if (g2 == null)
			throw new IllegalArgumentException("Object " + obj2 + " not contained in Percolator.");
		if (g1 != g2) {
			g1.addAll(g2);
			groups.remove(g2);
		}
	}
	
	public ArrayList<ArrayList<T>> getGroups() {
		return groups;
	}
}
