package kip.util


object RangeArray {
  /** A RangeArray that contains the results of some element computation a number of times. The range
   * step dx will be rounded up to be a multiple of xmax-min. */
  def fill[A: Manifest](xmin: Double, xmax: Double, dx: Double)(elem: => A): RangeArray[A] = {
    val n = ((xmax - xmin) / dx).toInt
    val bins = Array.fill[A](n)(elem)
    new RangeArray[A](xmin, xmax, bins)
  }
  
  def transpose[A: Manifest](c: Traversable[RangeArray[A]]): RangeArray[Array[A]] = {
    c.headOption match {
      case None => throw new IllegalArgumentException
      case Some(ra) => {
        val res = fill(ra.xmin, ra.xmax, ra.dx)(new collection.mutable.ArrayBuffer[A])
        for (r <- res.binCenters; a <- c) {
          res(r) += a(r)
        }
        res.map(_.toArray)
      }
    }
  }
}

class RangeArray[A](val xmin: Double, val xmax: Double, val elems: Array[A]) {
  val n = elems.size
  val dx = (xmax - xmin) / n
  
  def isDefinedAt(x: Double): Boolean = {
    xmin <= x && x <= xmax
  }
  
  def binCenterForIndex(i: Int): Double = {
    xmin + (i + 0.5) * dx
  }
  
  def binIndex(x: Double): Int = {
    if (!isDefinedAt(x))
      throw new IllegalArgumentException(
        "Argument %g out of range [%g, %g]" format (x, xmin, xmax))
    val i = ((x - xmin) / dx).toInt
    math.min(math.max(i, 0), n-1)
  }

  def binCenters[Double] = Array.tabulate(n)(binCenterForIndex _)
  
  def apply(x: Double): A = elems(binIndex(x))
  
  def update(x: Double, elem: A) {
    elems(binIndex(x)) = elem
  }
  
  def map[B: Manifest](f: A => B): RangeArray[B] = {
    new RangeArray(xmin, xmax, elems.map(f))
  }
  
  def exists(f: A => Boolean): Boolean = elems.exists(f)
  def foreach(f: A => Unit): Unit = elems.foreach(f)
}
