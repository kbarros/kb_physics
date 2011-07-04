package kip.math
package linalg

import com.sun.jna.ptr.IntByReference


trait DenseComplexMatrixImplicits {
  def c2r(a: Seq[Complex]): Array[Double] = DenseComplexMatrix.complexToDoubleArray(a)

  implicit def tuple2ToRow[T <% Complex](t: Tuple2[T, T]): DenseComplexMatrix = {
    new DenseComplexMatrix(1, 2, c2r(Seq(t._1, t._2)))
  }

  implicit def tuple3ToRow[T <% Complex](t: Tuple3[T, T, T]): DenseComplexMatrix = {
    new DenseComplexMatrix(1, 3, c2r(Seq(t._1, t._2, t._3)))
  }

  implicit def tuple4ToRow[T <% Complex](t: Tuple4[T, T, T, T]): DenseComplexMatrix = {
    new DenseComplexMatrix(1, 4, c2r(Seq(t._1, t._2, t._3, t._4)))
  }

  implicit def tuple5ToRow[T <% Complex](t: Tuple5[T, T, T, T, T]): DenseComplexMatrix = {
    new DenseComplexMatrix(1, 5, c2r(Seq(t._1, t._2, t._3, t._4, t._5)))
  }

  implicit def traversableToRow(t: Traversable[Complex]): DenseComplexMatrix = {
    new DenseComplexMatrix(1, t.size, c2r(t.toSeq))
  }

  implicit def arrayToRow(a: Array[Complex]): DenseComplexMatrix = {
    new DenseComplexMatrix(1, a.size, c2r(a))
  }
}

object DenseComplexMatrix extends DenseComplexMatrixImplicits {
  import Complex._
  
  def complexToDoubleArray(a: Seq[Complex]): Array[Double] = {
    val ret = new Array[Double](2*a.size)
    for (i <- a.indices) {
      ret(2*i+0) = a(i).re
      ret(2*i+1) = a(i).im
    }
    ret
  }

  def fill(numRows: Int, numCols: Int)(x: Complex) = {
    val elems = complexToDoubleArray(Array.fill(numRows*numCols)(x))
    new DenseComplexMatrix(numRows, numCols, elems)
  }
  
  def zeros(numRows: Int, numCols: Int) = {
    fill(numRows, numCols)(0)
  }
  
  def tabulate(numRows: Int, numCols: Int)(f: (Int, Int) => Complex) = {
    val m = zeros(numRows, numCols)
    for (j <- 0 until numCols; i <- 0 until numRows) m(i, j) = f(i, j)
    m
  }
  
  def eye(numRows: Int) = {
    val ret = zeros(numRows, numRows)
    for (i <- 0 until numRows)
      ret(i,i) = 1d
    ret
  }
  
  def column(elems: Complex*) = {
    new DenseComplexMatrix(elems.size, 1, complexToDoubleArray(elems))
  }
  
  def row(elems: Complex*) = {
    new DenseComplexMatrix(1, elems.size, complexToDoubleArray(elems))
  }
  
  def fromRows(row1: DenseComplexMatrix, rows: DenseComplexMatrix*) = {
    require(row1.numRows == 1 && rows.forall(_.numRows == 1))
    require(rows.forall(_.numCols == row1.numCols))
    
    val ret = zeros(1 + rows.size, row1.numCols)
    ret(0, ::) = row1
    for (i <- rows.indices)
      ret(i+1, ::) = rows(i)
    ret
  }
  
  // Native operations
  
  def blasZgemm(c: DenseComplexMatrix, a: DenseComplexMatrix, b: DenseComplexMatrix) {
    if (b.numCols == 1) {
      // TODO: optimize using cgemv
    }
    
    // c = alpha*a*b + beta*c
    Netlib.blas.zgemm("N", "N",
                      c.numRows, c.numCols, // dimension of return matrix
                      a.numCols, // dimension of summation index
                      Array(1.0, 0.0), // alpha 
                      a.data, a.numRows, // A matrix
                      b.data, b.numRows, // B matrix
                      Array(0.0, 0.0), // beta
                      c.data, c.numRows // C matrix
                    )
  }
  
  
  /** X := A \ V */
  def QRSolve(X : DenseComplexMatrix, A : DenseComplexMatrix, V : DenseComplexMatrix, transpose : Boolean) = {
    require(X.numRows == A.numCols, "Wrong number of rows in return value");
    require(X.numCols == V.numCols, "Wrong number of rows in return value");
    
    val nrhs = V.numCols;
    
    // allocate temporary solution matrix
    val Xtmp = DenseComplexMatrix.zeros(math.max(A.numRows, A.numCols), nrhs); // HOW COMPILE WITH REAL?
    val M = if (!transpose) A.numRows else A.numCols;
    for (j <- 0 until nrhs; i <- 0 until M) { Xtmp(i,j) = V(i,j): Complex; }

    val newData = A.data.clone();

    // query optimal workspace
    val queryWork = new Array[Double](2);
    val queryInfo = new IntByReference(0);
    Netlib.lapack.zgels(
      if (!transpose) "N" else "T",
      A.numRows, A.numCols, nrhs,
      newData, math.max(1,A.numRows),
      Xtmp.data, math.max(1,math.max(A.numRows,A.numCols)),
      queryWork, -1, queryInfo)
    
    // allocate workspace
    val workSize = 
      if (queryInfo.getValue != 0)
        math.max(1, math.min(A.numRows, A.numCols) + math.max(math.min(A.numRows, A.numCols), nrhs));
      else
        math.max(queryWork(0).toInt, 1);
    val work = new Array[Double](2*workSize); // multiply by two for re/im parts
    
    // compute factorization
    Netlib.lapack.zgels(
      if (!transpose) "N" else "T",
      A.numRows, A.numCols, nrhs,
      newData, math.max(1,A.numRows),
      Xtmp.data, math.max(1,math.max(A.numRows,A.numCols)),
      work, workSize, queryInfo);
    
    if (queryInfo.getValue< 0)
      throw new IllegalArgumentException;

    // extract solution
    val N = if (!transpose) A.numCols else A.numRows;
    for (j <- 0 until nrhs; i <- 0 until N) X(i,j) = Xtmp(i,j);

    X;
  }

}


class DenseComplexMatrix(val numRows: Int, val numCols: Int, val data: Array[Double]) {
  require(2*numRows*numCols == data.size)
  
  override def clone(): DenseComplexMatrix = {
    new DenseComplexMatrix(numRows, numCols, data.clone())
  }
  
  def checkKey(i: Int, j: Int) {
    require(0 <= i && i < numRows)
    require(0 <= j && j < numCols)
  }
  
  def index(i: Int, j: Int) = {
    checkKey(i, j)
    i + j*numRows // fortran column major convention
  }

  def apply(i: Int, j: Int): Complex = {
    val idx = index(i, j)
    Complex(data(2*idx+0), data(2*idx+1))
  }

  def apply(i: Int, _slice: Slice): DenseComplexMatrix = {
    DenseComplexMatrix.tabulate(1, numCols) { (_,j) => this(i,j) }
  }
  
  def apply(_slice: Slice, j: Int): DenseComplexMatrix = {
    DenseComplexMatrix.tabulate(numRows, 1) { (i,_) => this(i,j) }
  }

  def update(i: Int, j: Int, x: Complex) {
    val idx = index(i, j)
    data(2*idx+0) = x.re
    data(2*idx+1) = x.im
  }
  
  def update(i: Int, _slice: Slice, x: DenseComplexMatrix) {
    require(x.numRows == 1 && x.numCols == numCols)
    for (j <- 0 until numCols) this(i, j) = x(0, j)
  }

  def update(_slice: Slice, j: Int, x: DenseComplexMatrix) {
    require(x.numRows == numRows && x.numCols == 1)
    for (i <- 0 until numRows) this(i, j) = x(i, 0)
  }
  
  def transpose: DenseComplexMatrix = {
    DenseComplexMatrix.tabulate(numRows, numCols) { (i, j) => this(j, i) }
  }
  
  def conj : DenseComplexMatrix = {
    DenseComplexMatrix.tabulate(numRows, numCols) { (i, j) => this(i, j).conj }
  }
  
  def dagger : DenseComplexMatrix = {
    DenseComplexMatrix.tabulate(numRows, numCols) { (i, j) => this(j, i).conj }
  }
  
  def +(that: DenseComplexMatrix): DenseComplexMatrix = {
    require(numRows == that.numRows && numCols == that.numCols)
    DenseComplexMatrix.tabulate(numRows, numCols) { (i, j) => this(i, j) + that(i, j) }
  }

  def -(that: DenseComplexMatrix): DenseComplexMatrix = {
    require(numRows == that.numRows && numCols == that.numCols)
    DenseComplexMatrix.tabulate(numRows, numCols) { (i, j) => this(i, j) - that(i, j) }
  }

  def *(that: DenseComplexMatrix): DenseComplexMatrix = {
    require(numCols == that.numRows, "Size mismatch: cols %d != rows %d".format(numCols, that.numRows))
    val ret = DenseComplexMatrix.zeros(numRows, that.numCols)
    DenseComplexMatrix.blasZgemm(ret, this, that)
    // equivalent to:
    // for (i <- 0 until numRows;
    //      k <- 0 until numCols;
    //      j <- 0 until that.numCols) {
    //   ret(i, j) += this(i, k)*that(k, j)
    // }
    ret
  }
  
  def \(that: DenseComplexMatrix): DenseComplexMatrix = {
    require(numCols == that.numRows, "Size mismatch: cols %d != rows %d".format(numCols, that.numRows))
    require(that.numCols == 1)
    val ret = DenseComplexMatrix.zeros(numRows, 1)
    DenseComplexMatrix.QRSolve(ret, this, that, false)
    ret
  }

  def *(that: Complex): DenseComplexMatrix = {
    DenseComplexMatrix.tabulate(numRows, numCols) { (i, j) => this(i, j)*that }
  }

  def /(that: Double): DenseComplexMatrix = {
    DenseComplexMatrix.tabulate(numRows, numCols) { (i, j) => this(i, j)/that }
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

