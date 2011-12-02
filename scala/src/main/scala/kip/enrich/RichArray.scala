package kip.enrich

class RichDoubleArray(a: Array[Double]) {
  def +=[T <% Double](s: Seq[T]) {
    require(a.size == s.size)
    for (i <- 0 until a.size) {
      a(i) += s(i)
    }
  }
  
  def dot[T <% Double](s: Seq[T]) = {
    require(a.size == s.size)
    var ret = 0d
    for (i <- 0 until a.size) {
      ret += a(i) * s(i)
    }
    ret
  }
}


class RichFloatArray(a: Array[Float]) {
  def +=[T <% Float](s: Seq[T]) {
    require(a.size == s.size)
    for (i <- 0 until a.size) {
      a(i) += s(i)
    }
  }
  
  def dot[T <% Float](s: Seq[T]) = {
    require(a.size == s.size)
    var ret = 0d
    for (i <- 0 until a.size) {
      ret += a(i) * s(i)
    }
    ret
  }
}
