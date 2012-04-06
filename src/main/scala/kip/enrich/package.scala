package kip

import collection.TraversableLike
import java.io.File

package object enrich {
  implicit def enrichTraversableLike[Repr[a], A](s: Repr[A])(implicit ev: Repr[A] => TraversableLike[A, Repr[A]]) =
    new RichTraversableLike(s)
  
  implicit def enrichFile(f: File): RichFile = new RichFile(f)
  
  implicit def enrichDoubleArray(a: Array[Double]) = new RichDoubleArray(a)
  implicit def enrichFloatArray(a: Array[Float])   = new RichFloatArray(a)
}
