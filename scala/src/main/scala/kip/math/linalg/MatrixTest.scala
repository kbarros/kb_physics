package kip.math
package linalg

import com.sun.jna.ptr.IntByReference
import java.nio._

/*
 * FUTURE WORK:
 *
 * - Make ScalarData extend Scalar : Associate operations with type of Matrix
 * - LAPACK operations
 * - CanBuildFrom for matrices
 * - Sparse matrices
 * - Unify ops names
 *   - Generics for float precision
 * 
 * */

// ----------------------------------
// Scalar

trait Scalar[@specialized(Float, Double) A] {
  def add(a: A, b: A): A
  def sub(a: A, b: A): A
  def mul(a: A, b: A): A
  def div(a: A, b: A): A
  def negate(a: A): A
  def zero: A
}

trait DoubleScalar extends Scalar[Double] {
  def add(a: Double, b: Double): Double = a + b
  def sub(a: Double, b: Double): Double = a - b
  def mul(a: Double, b: Double): Double = a * b
  def div(a: Double, b: Double): Double = a / b
  def negate(a: Double): Double = -a
  def zero: Double = 0.0
}

trait ComplexScalar extends Scalar[Complex] {
  def add(a: Complex, b: Complex): Complex = a + b
  def sub(a: Complex, b: Complex): Complex = a - b
  def mul(a: Complex, b: Complex): Complex = a * b
  def div(a: Complex, b: Complex): Complex = a / b
  def negate(a: Complex): Complex = -a
  def zero: Complex = 0
}

/*
 * TODO: move into Scalar
object ScalarOps {
  trait ScalarOps[@specialized(Float, Double) A] {
    val lhs: A
    val n: Scalar[A]
    def +(rhs:A) = n.add(lhs, rhs)
    def -(rhs:A) = n.sub(lhs, rhs)
    def *(rhs:A) = n.mul(lhs, rhs)
    def /(rhs:A) = n.div(lhs, rhs)
    def unary_-() = n.negate(lhs)
  }
  implicit def infixScalarOps[@specialized(Float, Double) A: Scalar](a: A): ScalarOps[A] = new ScalarOps[A] {
    val lhs = a
    val n = implicitly[Scalar[A]]
  }
}
*/

object Scalar {
  // def scalar[@specialized(Float, Double) A: Scalar]() = implicitly[Scalar[A]]

  implicit object DoubleScalar extends DoubleScalar
  implicit object ComplexScalar extends ComplexScalar

  implicit object DoubleData extends DoubleData
//  implicit object ComplexData extends ComplexData
}


// ----------------------------------
// Data

trait Data[@specialized(Float, Double) A, B] extends Scalar[A] {
  def alloc(size: Int): B
  def copyTo(src: B, dst: B, start: Int, len: Int): Unit
  def copyElemTo(src: B, dst: B, i: Int): Unit
  def update(a: B, i: Int, x: A): Unit
  def apply(a: B, i: Int): A
  def size(a: B): Int
  def dispose(a: B): Unit
  
  def madd(a0: B, i0: Int, a1: B, i1: Int, a2: B, i2: Int) {
    val x0 = apply(a0, i0)
    val x1 = apply(a1, i1)
    val x2 = apply(a2, i2)
    update(a0, i0, add(x0, mul(x1, x2)))
  }
  
  def madd(a0: B, i0: Int, a1: B, i1: Int, x2: A) {
    val x0 = apply(a0, i0)
    val x1 = apply(a1, i1)
    update(a0, i0, add(x0, mul(x1, x2)))
  }
}

class DoubleData extends Data[Double, Array[Double]] with DoubleScalar {
  def alloc(size: Int) = new Array[Double](size)
  def copyTo(src: Array[Double], dst: Array[Double], start: Int, len: Int) { src.copyToArray(dst, start, len) }
  def copyElemTo(src: Array[Double], dst: Array[Double], i: Int) { dst(i) = src(i) }
  def update(a: Array[Double], i: Int, x: Double) { a(i) = x }
  def apply(a: Array[Double], i: Int): Double = a(i)
  def size(a: Array[Double]) = a.size
  def dispose(a: Array[Double]) { }
}

/*
class ComplexData(val data: Array[Double]) extends Data[Complex, ComplexData] {
  def scalar = Scalar.ComplexScalar
  def copyTo(that: ComplexData, start: Int, len: Int) { data.copyToArray(that.data, 2*start, 2*len) }
  def copyElemTo(that: ComplexData, i: Int) {
    that.data(2*i+0) = data(2*i+0)
    that.data(2*i+1) = data(2*i+1)
  }
  def update(i: Int, x: Complex) {
    data(2*i+0) = x.re
    data(2*i+1) = x.im
  }
  def apply(i: Int): Complex = Complex(data(2*i+0), data(2*i+1))
  def size = data.size/2
  def dispose() { }
  
  override def madd(i0: Int, a1: ComplexData, i1: Int, a2: ComplexData, i2: Int) {
    val a1_re = a1.data(2*i1+0)
    val a1_im = a1.data(2*i1+1)
    val a2_re = a2.data(2*i1+0)
    val a2_im = a2.data(2*i1+1)
    data(2*i0+0) += a1_re*a2_re - a1_im*a2_im
    data(2*i0+1) += a1_re*a2_im + a1_im*a2_re
  }
}
*/


// ----------------------------------
// Matrix impl

object MyDenseMatrix {
  def fill[@specialized(Float, Double) A, B, C <: Data[A, B]]
      (numRows: Int, numCols: Int)(x: A)
      (implicit ops: C): MyDenseMatrix[A, B, C] = {
    val data = ops.alloc(numRows*numCols)
    for (i <- 0 until ops.size(data)) ops.update(data, i, x)
    new MyDenseMatrix[A, B, C](numRows, numCols, data)
  }

  def zeros[@specialized(Float, Double) A, B, C <: Data[A, B]]
      (numRows: Int, numCols: Int)
      (implicit ops: C): MyDenseMatrix[A, B, C] = {
    fill(numRows, numCols)(ops.zero)
  }

  def tabulate[@specialized(Float, Double) A, B, C <: Data[A, B]]
      (numRows: Int, numCols: Int)(f: (Int, Int) => A)
      (implicit ops: C): MyDenseMatrix[A, B, C] = {
    val m = zeros[A, B, C](numRows, numCols)
    for (j <- 0 until numCols; i <- 0 until numRows) m(i, j) = f(i, j)
    m
  }
}


class MyDenseMatrix[@specialized(Float, Double) A, B, C <: Data[A, B]]
    (val numRows: Int, val numCols: Int, val data: B)
    (implicit ops: C) {

  require(numRows*numCols == ops.size(data))
  
  override def clone(): MyDenseMatrix[A, B, C] = {
    val newData: B = ops.alloc(ops.size(data))
    ops.copyTo(data, newData, 0, ops.size(data))
    new MyDenseMatrix(numRows, numCols, newData)
  }
  
  // def toComplex: DenseComplexMatrix = {
  //   DenseComplexMatrix.tabulate(numRows, numCols) { (i,j) => this(i,j) }
  // }
  
  def checkKey(i: Int, j: Int) {
    require(0 <= i && i < numRows)
    require(0 <= j && j < numCols)
  }
  
  def index(i: Int, j: Int) = {
    checkKey(i, j)
    i + j*numRows // fortran column major convention
  }

  def apply(i: Int, j: Int): A = {
    ops.apply(data, index(i, j))
  }
  
  def apply(i: Int, _slice: Slice): MyDenseMatrix[A, B, C] = {
    MyDenseMatrix.tabulate(1, numCols) { (_,j) => this(i,j) }
  }
  
  def apply(_slice: Slice, j: Int): MyDenseMatrix[A, B, C] = {
    MyDenseMatrix.tabulate(numRows, 1) { (i,_) => this(i,j) }
  }

  def update(i: Int, j: Int, x: A) {
    ops.update(data, index(i, j), x)
  }
  
  def update(i: Int, _slice: Slice, x: MyDenseMatrix[A, B, C]) {
    require(x.numRows == 1 && x.numCols == numCols)
    for (j <- 0 until numCols) this(i, j) = x(0, j)
  }

  def update(_slice: Slice, j: Int, x: MyDenseMatrix[A, B, C]) {
    require(x.numRows == numRows && x.numCols == 1)
    for (i <- 0 until numRows) this(i, j) = x(i, 0)
  }
  
  def transpose: MyDenseMatrix[A, B, C] = {
    MyDenseMatrix.tabulate(numRows, numCols) { (i, j) => this(j, i) }
  }
  
  def +(that: MyDenseMatrix[A, B, C]): MyDenseMatrix[A, B, C] = {
    require(numRows == that.numRows && numCols == that.numCols)
    MyDenseMatrix.tabulate(numRows, numCols) { (i, j) => ops.add(this(i, j), that(i, j)) }
  }

  def -(that: MyDenseMatrix[A, B, C]): MyDenseMatrix[A, B, C] = {
    require(numRows == that.numRows && numCols == that.numCols)
    MyDenseMatrix.tabulate(numRows, numCols) { (i, j) => ops.sub(this(i, j), that(i, j)) }
  }

  def *(that: MyDenseMatrix[A, B, C]): MyDenseMatrix[A, B, C] = {
    require(numCols == that.numRows, "Size mismatch: cols %d != rows %d".format(numCols, numRows))
    val ret = MyDenseMatrix.zeros[A, B, C](numRows, that.numCols)
//    MyDenseMatrix.blasDgemm(ret, this, that)
    for (i <- 0 until numRows;
         k <- 0 until numCols;
         j <- 0 until that.numCols) {
      // equivalent to: ret(i, j) += this(i, k)*that(k, j)
      ops.madd(ret.data, index(i, j), this.data, index(i, k), that.data, index(k, j))
    }
    ret
  }
  
  def \(that: MyDenseMatrix[A, B, C]): MyDenseMatrix[A, B, C] = {
    require(numCols == that.numRows, "Size mismatch: cols %d != rows %d".format(numCols, numRows))
    require(that.numCols == 1)
    val ret = MyDenseMatrix.zeros[A, B, C](numRows, 1)
//    MyDenseMatrix.QRSolve(ret, this, that, false)
    ret
  }

  def *(that: A): MyDenseMatrix[A, B, C] = {
    MyDenseMatrix.tabulate(numRows, numCols) { (i, j) => ops.mul(this(i, j), that) }
  }

  def /(that: A): MyDenseMatrix[A, B, C] = {
    MyDenseMatrix.tabulate(numRows, numCols) { (i, j) => ops.div(this(i, j), that) }
  }
  
  override def toString = {
    val sb = new StringBuilder()
    val elemWidth = 8
    def writeStr(str: String) {
      val spaces = Seq.fill(math.max(2, elemWidth-str.size))(' ').mkString
      sb.append(spaces)
      sb.append(str)
    }
    val maxRows = 10
    val maxCols = 10
    for (i <- 0 until math.min(numRows, maxRows)) {
      for (j <- 0 until math.min(numCols, maxCols)) {
        writeStr(this(i, j).toString)
      }
      if (i == 0 && numCols > maxCols)
        sb.append(" ... (%d Cols)".format(numCols))
      sb.append("\n")
    }
    if (numRows > maxRows)
      sb.append(" ... (%d Rows)".format(numRows))
    sb.toString
  }
}
