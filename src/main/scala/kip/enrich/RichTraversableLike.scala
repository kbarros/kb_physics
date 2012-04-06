package kip.enrich

import collection.TraversableLike
import collection.generic.CanBuildFrom


class RichTraversableLike[Repr[a], A](s: Repr[A])(implicit ev: Repr[A] => TraversableLike[A, Repr[A]]) {
  val ths = ev(s)

  def foldMapLeft[B, That](z: B)(op: (B, A) => B)(implicit bf: CanBuildFrom[Repr[A], B, That]): That = {
    val b = bf(ths.repr)
    b.sizeHint(ths)
    var _z = z
    for (x <- ths) {
      _z = op(_z, x)
      b += _z
    }
    b.result
  }
}
