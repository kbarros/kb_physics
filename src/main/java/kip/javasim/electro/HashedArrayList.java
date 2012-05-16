package kip.javasim.electro;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

public class HashedArrayList<T> implements Iterable<T> {
	private ArrayList<T> list;
	private HashMap<T, Integer> hash;
	
	public HashedArrayList() {
		list = new ArrayList<T>();
		hash = new HashMap<T,Integer>();
	}
	
	public void add(T val) {
		if (hash.containsKey(val)) {
			System.err.println("List already contains " + val);
		}
		else {
			hash.put(val, list.size());
			list.add(val);
		}
	}
	
	public void replace(T before, T after) {
		int index = hash.get(before);
		if (list.get(index) != before)
			throw new IllegalStateException();
		hash.remove(before);
		hash.put(after, index);
		list.set(index, after);
	}
	
	public int size() {
		return list.size();
	}
	
	public T get(int index) {
		return list.get(index);
	}
	
	public Iterator<T> iterator() {
		return list.iterator();
	}
}
