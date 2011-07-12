package kip.math.linalg3

import kip.math.Complex
import kip.math.linalg2.{ScalarOps, ScalarData}

trait Matrix[A, Raw, +Repr[C, D] <: Matrix[C, D, Repr]] { self: Repr[A, Raw] =>
  val scalar: ScalarOps[A]
  
  def numRows: Int
  def numCols: Int
  def apply(i: Int, j: Int): A
  def update(i: Int, j: Int, x: A)
  def indices: Iterator[(Int, Int)]

  def clone[That[C,D] >: Repr[C,D]](implicit mb: MatrixBuilder[A, Raw, That]): That[A, Raw] = {
    val ret = mb.zeros(numRows, numCols)
//    copyTo(ret)
    ret
  }
  
  def +[Repr2[C, D] <: Matrix[C, D, Repr2], Repr3[C, D] <: Matrix[C, D, Repr3]]
        (that: Repr2[A, Raw])
        (implicit ma: MatrixAdder[A, Raw, Repr, Repr2, Repr3],
                  mb: MatrixBuilder[A, Raw, Repr3]): Repr3[A, Raw] = {
    require(numCols == that.numRows)
    val ret = mb.zeros(numRows, numCols)
    ma.addInPlace(false, this, false, that, ret)
    ret
  }

  def get: Repr[A, Raw] = this

}

object MatrixAdder extends DenseAdders

trait MatrixAdder[A, Raw, -Repr1[_, _], -Repr2[_, _], -Repr3[_, _]] {
  def addInPlace(neg1: Boolean, m1: Repr1[A, Raw], neg2: Boolean, m2: Repr2[A, Raw], ret: Repr3[A, Raw])
}

object MatrixBuilder extends DenseBuilders

trait MatrixBuilder[A, Raw, +Repr[_, _]] {
  def zeros(numRows: Int, numCols: Int): Repr[A, Raw]
}


object Test extends App {
  val m1 = MatrixBuilder.dense[Float, Float].zeros(4, 4)
  val m2 = MatrixBuilder.denseRow[Float, Float].zeros(4, 4)
  val m3 = m1+m2
  println(m3)

//  println(implicitly[MatrixAdder[Float, Float, DenseRow, DenseRow, Dense]])
}




/*
// why compiles?
trait Matrix[-Repr <: Matrix[Repr]]
trait Dense extends Matrix[Dense]
trait DenseRow extends Matrix[DenseRow] with Dense
*/
