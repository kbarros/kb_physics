package kip.math
package linalg

import com.sun.jna.ptr.IntByReference
import java.nio._

/*
 * FUTURE WORK:
 *
 * - LAPACK operations
 * - CanBuildFrom for matrices
 * - Sparse matrices
 * 
 * */

object DenseMatrix {
  implicit object RealDbl    extends Builder[Double,  Array[Double]] ()(ScalarData.RealDblData)
  implicit object RealFlt    extends Builder[Float,   Array[Float]]  ()(ScalarData.RealFltData)
  implicit object ComplexDbl extends Builder[Complex, Array[Double]] ()(ScalarData.ComplexDblData)
  implicit object ComplexFlt extends Builder[Complex, Array[Float]]  ()(ScalarData.ComplexFltData)

  class Builder[@specialized(Float, Double) A, B](implicit scalar: ScalarData[A, B]) {
    def fill(numRows: Int, numCols: Int)(x: A): DenseMatrix[A, B] = {
      val data = scalar.alloc(numRows*numCols)
      for (i <- 0 until scalar.size(data)) scalar.update(data, i, x)
      new DenseMatrix[A, B](numRows, numCols, data)
    }

    def zeros(numRows: Int, numCols: Int): DenseMatrix[A, B] = {
      fill(numRows, numCols)(scalar.zero)
    }

    def tabulate(numRows: Int, numCols: Int)(f: (Int, Int) => A): DenseMatrix[A, B] = {
      val m = zeros(numRows, numCols)
      for (j <- 0 until numCols; i <- 0 until numRows) m(i, j) = f(i, j)
      m
    }

    def eye(numRows: Int): DenseMatrix[A, B] = {
      tabulate(numRows, numRows) { (i, j) => if (i == j) scalar.one else scalar.zero }
    }
    
    def column(elems: A*): DenseMatrix[A, B] = {
      val data = scalar.alloc(elems.size)
      for (i <- 0 until elems.size) scalar.update(data, i, elems(i))
      new DenseMatrix(elems.size, 1, data)
    }

    def row(elems: A*): DenseMatrix[A, B] = {
      val data = scalar.alloc(elems.size)
      for (i <- 0 until elems.size) scalar.update(data, i, elems(i))
      new DenseMatrix(1, elems.size, data)
    }

    def fromRows(row1: DenseMatrix[A, B], rows: DenseMatrix[A, B]*): DenseMatrix[A, B] = {
      require(row1.numRows == 1 && rows.forall(_.numRows == 1))
      require(rows.forall(_.numCols == row1.numCols))
      val ret = zeros(1 + rows.size, row1.numCols)
      ret(0, ::) = row1
      for (i <- rows.indices)
        ret(i+1, ::) = rows(i): DenseMatrix[A, B]
      ret
    }
  }
}


class DenseMatrix[@specialized(Float, Double) A, B]
    (val numRows: Int, val numCols: Int, val data: B)
    (implicit scalar: ScalarData[A, B]) {

  require(numRows*numCols == scalar.size(data))
  val matrix = new DenseMatrix.Builder
  
  override def clone(): DenseMatrix[A, B] = {
    val newData: B = scalar.alloc(scalar.size(data))
    scalar.copyTo(data, newData, 0, scalar.size(data))
    new DenseMatrix(numRows, numCols, newData)
  }
  
  // TODO: Introduce "Repr" parameter to get most specific return type
  def map[A0, B0](f: A => A0)(implicit builder: DenseMatrix.Builder[A0, B0]): DenseMatrix[A0, B0] = {
    builder.tabulate(numRows, numCols) { (i, j) => f(this(i, j)) }
  }
  
  def checkKey(i: Int, j: Int) {
    require(0 <= i && i < numRows)
    require(0 <= j && j < numCols)
  }
  
  def index(i: Int, j: Int) = {
    println("getting key %d %d for dim %d %d\n".format(i, j, numRows, numCols))
    checkKey(i, j)
    i + j*numRows // fortran column major convention
  }

  def apply(i: Int, j: Int): A = {
    scalar.apply(data, index(i, j))
  }
  
  def apply(i: Int, _slice: Slice): DenseMatrix[A, B] = {
    matrix.tabulate(1, numCols) { (_,j) => this(i,j) }
  }
  
  def apply(_slice: Slice, j: Int): DenseMatrix[A, B] = {
    matrix.tabulate(numRows, 1) { (i,_) => this(i,j) }
  }

  def update(i: Int, j: Int, x: A) {
    scalar.update(data, index(i, j), x)
  }
  
  def update(i: Int, _slice: Slice, x: DenseMatrix[A, B]) {
    require(x.numRows == 1 && x.numCols == numCols)
    for (j <- 0 until numCols) this(i, j) = x(0, j)
  }

  def update(_slice: Slice, j: Int, x: DenseMatrix[A, B]) {
    require(x.numRows == numRows && x.numCols == 1)
    for (i <- 0 until numRows) this(i, j) = x(i, 0)
  }
  
  def +(that: DenseMatrix[A, B]): DenseMatrix[A, B] = {
    require(numRows == that.numRows && numCols == that.numCols)
    matrix.tabulate(numRows, numCols) { (i, j) => scalar.add(this(i, j), that(i, j)) }
  }

  def -(that: DenseMatrix[A, B]): DenseMatrix[A, B] = {
    require(numRows == that.numRows && numCols == that.numCols)
    matrix.tabulate(numRows, numCols) { (i, j) => scalar.sub(this(i, j), that(i, j)) }
  }

  def *(that: DenseMatrix[A, B]): DenseMatrix[A, B] = {
    require(numCols == that.numRows, "Size mismatch: cols %d != rows %d".format(numCols, numRows))
    val ret = matrix.zeros(numRows, that.numCols)
//    matrix.blasDgemm(ret, this, that)
    for (i <- 0 until numRows;
         k <- 0 until numCols;
         j <- 0 until that.numCols) {
      // equivalent to: ret(i, j) += this(i, k)*that(k, j)
      scalar.madd(ret.data, ret.index(i, j), data, index(i, k), that.data, that.index(k, j))
    }
    ret
  }
  
  def \(that: DenseMatrix[A, B]): DenseMatrix[A, B] = {
    require(numCols == that.numRows, "Size mismatch: cols %d != rows %d".format(numCols, numRows))
    require(that.numCols == 1)
    val ret = matrix.zeros(numRows, 1)
//    matrix.QRSolve(ret, this, that, false)
    ret
  }

  def *(that: A): DenseMatrix[A, B] = {
    matrix.tabulate(numRows, numCols) { (i, j) => scalar.mul(this(i, j), that) }
  }

  def /(that: A): DenseMatrix[A, B] = {
    matrix.tabulate(numRows, numCols) { (i, j) => scalar.div(this(i, j), that) }
  }

  def unary_- : DenseMatrix[A, B] = {
    matrix.tabulate(numRows, numCols) { (i, j) => scalar.neg(this(i, j)) }
  }
  
  def tran: DenseMatrix[A, B] = {
    matrix.tabulate(numCols, numRows) { (i, j) => this(j, i) }
  }
  
  def conj: DenseMatrix[A, B] = {
    matrix.tabulate(numRows, numCols) { (i, j) => scalar.conj(this(i, j)) }
  }
  
  def dag: DenseMatrix[A, B] = {
    matrix.tabulate(numCols, numRows) { (i, j) => scalar.conj(this(j, i)) }
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

