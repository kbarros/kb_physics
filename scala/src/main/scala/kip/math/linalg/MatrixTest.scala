package kip.math
package linalg

import com.sun.jna.ptr.IntByReference
import java.nio._

/*
 * FUTURE WORK:
 *
 * - Nicer syntax, type inference
 * - LAPACK operations
 * - CanBuildFrom for matrices
 * - Sparse matrices
 * 
 * */

// ----------------------------------
// Matrix impl

object MyDenseMatrix {
  def fill[@specialized(Float, Double) A, B, C <: ScalarData[A, B]]
      (numRows: Int, numCols: Int)(x: A)
      (implicit ops: C): MyDenseMatrix[A, B, C] = {
    val data = ops.alloc(numRows*numCols)
    for (i <- 0 until ops.size(data)) ops.update(data, i, x)
    new MyDenseMatrix[A, B, C](numRows, numCols, data)
  }

  def zeros[@specialized(Float, Double) A, B, C <: ScalarData[A, B]]
      (numRows: Int, numCols: Int)
      (implicit ops: C): MyDenseMatrix[A, B, C] = {
    fill(numRows, numCols)(ops.zero)
  }

  def tabulate[@specialized(Float, Double) A, B, C <: ScalarData[A, B]]
      (numRows: Int, numCols: Int)(f: (Int, Int) => A)
      (implicit ops: C): MyDenseMatrix[A, B, C] = {
    val m = zeros[A, B, C](numRows, numCols)
    for (j <- 0 until numCols; i <- 0 until numRows) m(i, j) = f(i, j)
    m
  }
}


class MyDenseMatrix[@specialized(Float, Double) A, B, C <: ScalarData[A, B]]
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


/*
trait ScalarTypeclass {
  type T
  def foo(x: T): Unit
}
object Test {
  def test[C <: ScalarTypeclass](scalar: C, x: C#T) = scalar.foo(x)
}

*/

