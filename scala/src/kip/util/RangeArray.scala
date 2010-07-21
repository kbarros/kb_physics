package kip.util


object RangeArray {
  def fill[A: Manifest](xmin: Double, xmax: Double, dx: Double)(elem: => A): RangeArray[A] = {
    val n = ((xmax - xmin) / dx).toInt
    val bins = Array.fill[A](n)(elem)
    new RangeArray[A](xmin, xmax, bins)
  }
  
  def transpose[A, T[_] <: Traversable[_]](c: T[RangeArray[A]]): RangeArray[T[A]] = {
    // TODO: figure out how to implement this :-)
    null: RangeArray[T[A]]
  }
}

class RangeArray[A](val xmin: Double, val xmax: Double, val elems: Array[A]) {
  val n = elems.size
  val dx = (xmax - xmin) / n
  
  
  def isDefinedAt(x: Double): Boolean = {
    xmin <= x && x <= xmax
  }
  
  def positionOfIndex(i: Int): Double = {
    xmin + (i + 0.5) * dx
  }
  
  def indexOfPosition(x: Double): Int = {
    if (!isDefinedAt(x))
      throw new IllegalArgumentException(
        "Argument %g out of range [%g, %g]" format (x, xmin, xmax))
    val i = ((x - xmin) / dx).toInt
    math.min(math.max(i, 0), n-1)
  }
  
  def positions: Array[Double] = Array.tabulate(n)(positionOfIndex _)
  
  def apply(x: Double): A = {
    elems(indexOfPosition(x))
  }
  
  def update(x: Double, elem: A) {
    elems(indexOfPosition(x)) = elem
  }

  def map[B: Manifest](f: A => B): RangeArray[B] = {
    new RangeArray(xmin, xmax, elems.map(f))
  }
}
