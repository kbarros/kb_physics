package kip.math
package linalg

import com.sun.jna.ptr.IntByReference


trait DenseMatrixImplicits {
  implicit def tuple2ToRow[T <% Double](t: Tuple2[T, T]): DenseMatrix = {
    new DenseMatrix(1, 2, Array(t._1, t._2))
  }
  
  implicit def tuple3ToRow[T <% Double](t: Tuple3[T, T, T]): DenseMatrix = {
    new DenseMatrix(1, 3, Array(t._1, t._2, t._3))
  }

  implicit def tuple4ToRow[T <% Double](t: Tuple4[T, T, T, T]): DenseMatrix = {
    new DenseMatrix(1, 4, Array(t._1, t._2, t._3, t._4))
  }

  implicit def tuple5ToRow[T <% Double](t: Tuple5[T, T, T, T, T]): DenseMatrix = {
    new DenseMatrix(1, 5, Array(t._1, t._2, t._3, t._4, t._5))
  }

  implicit def traversableToRow(t: Traversable[Double]): DenseMatrix = {
    new DenseMatrix(1, t.size, t.toArray)
  }

  implicit def arrayToRow(a: Array[Double]): DenseMatrix = {
    new DenseMatrix(1, a.size, a)
  }
}


object DenseMatrix extends DenseMatrixImplicits {

  def fill(numRows: Int, numCols: Int)(x: Double) = {
    val elems = Array.fill(numRows*numCols)(x)
    new DenseMatrix(numRows, numCols, elems)
  }
  
  def zeros(numRows: Int, numCols: Int) = {
    fill(numRows, numCols)(0)
  }
  
  def tabulate(numRows: Int, numCols: Int)(f: (Int, Int) => Double) = {
    val m = zeros(numRows, numCols)
    for (j <- 0 until numCols; i <- 0 until numRows) m(i, j) = f(i, j)
    m
  }

  def eye(numRows: Int) = {
    val ret = zeros(numRows, numRows)
    for (i <- 0 until numRows)
      ret(i,i) = 1
    ret
  }
  
  def column(elems: Double*) = {
    new DenseMatrix(elems.size, 1, elems.toArray)
  }
    
  def row(elems: Double*) = {
    new DenseMatrix(1, elems.size, elems.toArray)
  }
  
  def fromRows(row1: DenseMatrix, rows: DenseMatrix*) = {
    require(row1.numRows == 1 && rows.forall(_.numRows == 1))
    require(rows.forall(_.numCols == row1.numCols))
    
    val ret = zeros(1 + rows.size, row1.numCols)
    ret(0, ::) = row1
    for (i <- rows.indices)
      ret(i+1, ::) = rows(i)
    ret
  }
  
  // Native operations
  
  // c = alpha*a*b + beta*c
  def blasDgemm(c: DenseMatrix, a: DenseMatrix, b: DenseMatrix) {
    Netlib.blas.dgemm("N", "N",
                      c.numRows, c.numCols, // dimension of return matrix
                      a.numCols, // dimension of summation index
                      1.0, // alpha 
                      a.data, a.numRows, // A matrix
                      b.data, b.numRows, // B matrix
                      0.0, // beta
                      c.data, c.numRows // C matrix
                    )
  }
  

  /** X := A \ V */
  def QRSolve(X : DenseMatrix, A : DenseMatrix, V : DenseMatrix, transpose : Boolean) = {
    require(X.numRows == A.numCols, "Wrong number of rows in return value");
    require(X.numCols == V.numCols, "Wrong number of rows in return value");

    val nrhs = V.numCols;
    
    // allocate temporary solution matrix
    val Xtmp = DenseMatrix.zeros(math.max(A.numRows, A.numCols), nrhs);
    val M = if (!transpose) A.numRows else A.numCols;
    for (j <- 0 until nrhs; i <- 0 until M) { Xtmp(i,j) = V(i,j); }

    val newData = A.data.clone();

    // query optimal workspace
    val queryWork = new Array[Double](1);
    val queryInfo = new IntByReference(0);
    Netlib.lapack.dgels(
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
    val work = new Array[Double](workSize);

    // compute factorization
    Netlib.lapack.dgels(
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

  def eig(m : DenseMatrix): (DenseMatrix, DenseMatrix, DenseMatrix) = {
    require(m.numRows == m.numCols, "Operation 'eig' requires square matrix")

    val n = m.numRows;

    // Allocate space for the decomposition
    var Wr = DenseMatrix.zeros(n, 1);
    var Wi = DenseMatrix.zeros(n, 1);
    var Vr = DenseMatrix.zeros(n, n);

    // Find the needed workspace
    val worksize = Array.ofDim[Double](1);
    val info = new IntByReference(0);
    Netlib.lapack.dgeev(
      "N", "V", n,
      Array.empty[Double], math.max(1,n),
      Array.empty[Double], Array.empty[Double],
      Array.empty[Double], math.max(1,n),
      Array.empty[Double], math.max(1,n),
      worksize, -1, info);

    // Allocate the workspace
    val lwork: Int = if (info.getValue != 0)
      math.max(1,4*n);
    else
      math.max(1,worksize(0).toInt);
    val work = Array.ofDim[Double](lwork);

    // Factor it!
    Netlib.lapack.dgeev(
      "N", "V", n,
      m.clone().data, math.max(1,n),
      Wr.data, Wi.data,
      Array.empty[Double], math.max(1,n),
      Vr.data, math.max(1,n),
      work, lwork, info);

    require(info.getValue >= 0, "Error in dgeev argument %d".format(-info.getValue))
    require(info.getValue <= 0, "Not converged dgeev; only %d of %d eigenvalues computed".format(info.getValue, m.numRows))
    
    (Wr, Wi, Vr)
  }
}


class DenseMatrix(val numRows: Int, val numCols: Int, val data: Array[Double]) extends Cloneable {
  require(numRows*numCols == data.size)
  
  override def clone(): DenseMatrix = {
    new DenseMatrix(numRows, numCols, data.clone())
  }
  
  def toComplex: DenseComplexMatrix = {
    DenseComplexMatrix.tabulate(numRows, numCols) { (i,j) => this(i,j) }
  }
  
  def checkKey(i: Int, j: Int) {
    require(0 <= i && i < numRows)
    require(0 <= j && j < numCols)
  }
  
  def index(i: Int, j: Int) = {
    checkKey(i, j)
    i + j*numRows // fortran column major convention
  }

  def apply(i: Int, j: Int): Double = {
    data(index(i, j))
  }
  
  def apply(i: Int, _slice: Slice): DenseMatrix = {
    DenseMatrix.tabulate(1, numCols) { (_,j) => this(i,j) }
  }
  
  def apply(_slice: Slice, j: Int): DenseMatrix = {
    DenseMatrix.tabulate(numRows, 1) { (i,_) => this(i,j) }
  }

  def update(i: Int, j: Int, x: Double) {
    data(index(i, j)) = x
  }
  
  def update(i: Int, _slice: Slice, x: DenseMatrix) {
    require(x.numRows == 1 && x.numCols == numCols)
    for (j <- 0 until numCols) this(i, j) = x(0, j)
  }

  def update(_slice: Slice, j: Int, x: DenseMatrix) {
    require(x.numRows == numRows && x.numCols == 1)
    for (i <- 0 until numRows) this(i, j) = x(i, 0)
  }
  
  def transpose: DenseMatrix = {
    DenseMatrix.tabulate(numRows, numCols) { (i, j) => this(j, i) }
  }
  
  def +(that: DenseMatrix): DenseMatrix = {
    require(numRows == that.numRows && numCols == that.numCols)
    DenseMatrix.tabulate(numRows, numCols) { (i, j) => this(i, j) + that(i, j) }
  }

  def -(that: DenseMatrix): DenseMatrix = {
    require(numRows == that.numRows && numCols == that.numCols)
    DenseMatrix.tabulate(numRows, numCols) { (i, j) => this(i, j) - that(i, j) }
  }

  def *(that: DenseMatrix): DenseMatrix = {
    require(numCols == that.numRows, "Size mismatch: cols %d != rows %d".format(numCols, numRows))
    val ret = DenseMatrix.zeros(numRows, that.numCols)
    DenseMatrix.blasDgemm(ret, this, that)
    // equivalent to:
    // for (i <- 0 until numRows;
    //      k <- 0 until numCols;
    //      j <- 0 until that.numCols) {
    //   ret(i, j) += this(i, k)*that(k, j)
    // }
    ret
  }
  
  def \(that: DenseMatrix): DenseMatrix = {
    require(numCols == that.numRows, "Size mismatch: cols %d != rows %d".format(numCols, numRows))
    require(that.numCols == 1)
    val ret = DenseMatrix.zeros(numRows, 1)
    DenseMatrix.QRSolve(ret, this, that, false)
    ret
  }

  def *(that: Double): DenseMatrix = {
    DenseMatrix.tabulate(numRows, numCols) { (i, j) => this(i, j)*that }
  }

  def /(that: Double): DenseMatrix = {
    DenseMatrix.tabulate(numRows, numCols) { (i, j) => this(i, j)/that }
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

