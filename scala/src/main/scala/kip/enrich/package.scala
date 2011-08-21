package kip

import collection.TraversableLike

package object enrich {
  implicit def enrichTraversableLike[Repr[a], A](s: Repr[A])(implicit ev: Repr[A] => TraversableLike[A, Repr[A]]) =
    new RichTraversableLike(s)
}
