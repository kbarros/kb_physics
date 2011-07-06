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
 * - Generics for float precision
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
  // def factory[B: Factory]() = implicitly[Factory[B]]

  implicit object DoubleScalar extends DoubleScalar
  implicit object DoubleDataFactory extends DoubleDataFactory

  implicit object ComplexScalar extends ComplexScalar
  implicit object ComplexDataFactory extends ComplexDataFactory
}


// ----------------------------------
// Data

trait Data[@specialized(Float, Double) A, B <: Data[A, B]] {
  implicit def scalar: Scalar[A]
  def data: Any
  def copyTo(that: B, start: Int, len: Int): Unit
  def copyElemTo(that: B, i: Int): Unit
  def update(i: Int, x: A): Unit
  def apply(i: Int): A
  def size: Int
  def dispose(): Unit
  
  def madd(i0: Int, a1: B, i1: Int, a2: B, i2: Int) {
    this(i0) = scalar.add(this(i0), scalar.mul(a1(i1), a2(i2)))
  }
  
  def madd(i0: Int, a1: B, i1: Int, a2: A) {
    this(i0) = scalar.add(this(i0), scalar.mul(a1(i1), a2))
  }
}

class DoubleData(val data: Array[Double]) extends Data[Double, DoubleData] {
  def scalar = Scalar.DoubleScalar
  def copyTo(that: DoubleData, start: Int, len: Int) { data.copyToArray(that.data, start, len) }
  def copyElemTo(that: DoubleData, i: Int) { that.data(i) = data(i) }
  def update(i: Int, x: Double) { data(i) = x }
  def apply(i: Int): Double = data(i)
  def size = data.size
  def dispose() { }
}

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


// ----------------------------------
// DataFactory

trait DataFactory[B] {
  def alloc(size: Int): B
}

trait DoubleDataFactory extends DataFactory[DoubleData] {
  def alloc(size: Int): DoubleData = {
    new DoubleData(new Array[Double](size))
  }
}

trait ComplexDataFactory extends DataFactory[ComplexData] {
  def alloc(size: Int): ComplexData = {
    new ComplexData(new Array[Double](2*size))
  }
}


// ----------------------------------
// Matrix impl


object MyDenseMatrix {
  def fill[@specialized(Float, Double) A, B <: Data[A, B]: DataFactory]
      (numRows: Int, numCols: Int)(x: A): MyDenseMatrix[A, B] = {
    val data = implicitly[DataFactory[B]].alloc(numRows*numCols)
    for (i <- 0 until data.size) data(i) = x
    new MyDenseMatrix[A, B](numRows, numCols, data)
  }

  def zeros[@specialized(Float, Double) A, B <: Data[A, B]: DataFactory]
      (numRows: Int, numCols: Int): MyDenseMatrix[A, B] = {
    val data = implicitly[DataFactory[B]].alloc(numRows*numCols)
    for (i <- 0 until data.size) data(i) = data.scalar.zero
    new MyDenseMatrix[A, B](numRows, numCols, data)
  }

  def tabulate[@specialized(Float, Double) A, B <: Data[A, B]: DataFactory]
      (numRows: Int, numCols: Int) (f: (Int, Int) => A): MyDenseMatrix[A, B] = {
    val m = zeros[A, B](numRows, numCols)
    for (j <- 0 until numCols; i <- 0 until numRows) m(i, j) = f(i, j)
    m
  }
}

class MyDenseMatrix[@specialized(Float, Double) A, B <: Data[A, B]: DataFactory]
    (val numRows: Int, val numCols: Int, val data: B) {

  def factory = implicitly[DataFactory[B]]
  
  require(numRows*numCols == data.size)
  
  override def clone(): MyDenseMatrix[A, B] = {
    val newData: B = factory.alloc(data.size)
    data.copyTo(newData, 0, data.size)
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
    data(index(i, j))
  }
  
  def apply(i: Int, _slice: Slice): MyDenseMatrix[A, B] = {
    MyDenseMatrix.tabulate(1, numCols) { (_,j) => this(i,j) }
  }
  
  def apply(_slice: Slice, j: Int): MyDenseMatrix[A, B] = {
    MyDenseMatrix.tabulate(numRows, 1) { (i,_) => this(i,j) }
  }

  def update(i: Int, j: Int, x: A) {
    data(index(i, j)) = x
  }
  
  def update(i: Int, _slice: Slice, x: MyDenseMatrix[A, B]) {
    require(x.numRows == 1 && x.numCols == numCols)
    for (j <- 0 until numCols) this(i, j) = x(0, j)
  }

  def update(_slice: Slice, j: Int, x: MyDenseMatrix[A, B]) {
    require(x.numRows == numRows && x.numCols == 1)
    for (i <- 0 until numRows) this(i, j) = x(i, 0)
  }
  
  def transpose: MyDenseMatrix[A, B] = {
    MyDenseMatrix.tabulate(numRows, numCols) { (i, j) => this(j, i) }
  }
  
  def +(that: MyDenseMatrix[A, B]): MyDenseMatrix[A, B] = {
    require(numRows == that.numRows && numCols == that.numCols)
    MyDenseMatrix.tabulate(numRows, numCols) { (i, j) => data.scalar.add(this(i, j), that(i, j)) }
  }

  def -(that: MyDenseMatrix[A, B]): MyDenseMatrix[A, B] = {
    require(numRows == that.numRows && numCols == that.numCols)
    MyDenseMatrix.tabulate(numRows, numCols) { (i, j) => data.scalar.sub(this(i, j), that(i, j)) }
  }

  def *(that: MyDenseMatrix[A, B]): MyDenseMatrix[A, B] = {
    require(numCols == that.numRows, "Size mismatch: cols %d != rows %d".format(numCols, numRows))
    val ret = MyDenseMatrix.zeros[A, B](numRows, that.numCols)
//    MyDenseMatrix.blasDgemm(ret, this, that)
    for (i <- 0 until numRows;
         k <- 0 until numCols;
         j <- 0 until that.numCols) {
      // equivalent to: ret(i, j) += this(i, k)*that(k, j)
      ret.data.madd(index(i, j), this.data, index(i, k), that.data, index(k, j))
    }
    ret
  }
  
  def \(that: MyDenseMatrix[A, B]): MyDenseMatrix[A, B] = {
    require(numCols == that.numRows, "Size mismatch: cols %d != rows %d".format(numCols, numRows))
    require(that.numCols == 1)
    val ret = MyDenseMatrix.zeros[A, B](numRows, 1)
//    MyDenseMatrix.QRSolve(ret, this, that, false)
    ret
  }

  def *(that: A): MyDenseMatrix[A, B] = {
    MyDenseMatrix.tabulate(numRows, numCols) { (i, j) => data.scalar.mul(this(i, j), that) }
  }

  def /(that: A): MyDenseMatrix[A, B] = {
    MyDenseMatrix.tabulate(numRows, numCols) { (i, j) => data.scalar.div(this(i, j), that) }
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
