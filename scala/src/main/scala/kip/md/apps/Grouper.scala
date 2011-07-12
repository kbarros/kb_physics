package kip.md.apps

import scala.collection.mutable


object DisjointSet {
  def findRoot[A](x: DisjointSet[A]): DisjointSet[A] = {
    if (x.parent == x)
      x
    else {
      x.parent = findRoot(x.parent)
      x.parent
    }
  }

  def union[A](x: DisjointSet[A], y: DisjointSet[A]): Unit = {
    val xRoot = findRoot(x)
    val yRoot = findRoot(y)
    if (xRoot != yRoot) {
      // x and y are not already in same set. Merge them.
      if (xRoot.rank < yRoot.rank)
        xRoot.parent = yRoot
      else if (xRoot.rank > yRoot.rank)
        yRoot.parent = xRoot
      else {
        yRoot.parent = xRoot
        xRoot.rank = xRoot.rank + 1
      }
    }
  }
}

class DisjointSet[A](val elem: A) {
  var parent: DisjointSet[A] = this
  var rank: Int = 0
}


class Grouper[A](elems: Seq[A]) {
  val elemToNode = (elems zip elems.map(new DisjointSet(_))).toMap
  def bond(x: A, y: A): Unit = DisjointSet.union(elemToNode(x), elemToNode(y))
  def groupRep(x: A): DisjointSet[A] = DisjointSet.findRoot(elemToNode(x))
  def groupMap: Map[DisjointSet[A], Seq[A]] = elems.groupBy(x => DisjointSet.findRoot(elemToNode(x)))
}
