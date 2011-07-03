package kip.math.linalg

import kip.netlib.{BlasLibrary => Blas}
import kip.netlib.BlasLibrary.{INSTANCE => blas}


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
    val elems = Array.fill(numRows*numCols)(0d)
    new DenseMatrix(numRows, numCols, elems)
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
  
  def blasDgemm(c: DenseMatrix, a: DenseMatrix, b: DenseMatrix) {
    if (b.numCols == 1) {
      // TODO: optimize using dgemv
    }
    
    // c = alpha*a*b + beta*c
    blas.cblas_dgemm(Blas.CblasColMajor,
                     Blas.CblasNoTrans, Blas.CblasNoTrans,
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

    import com.sun.jna.ptr.IntByReference
    implicit def intToIntByReference(a: Int) = new IntByReference(a)
    
    val nrhs = V.numCols;
    
    // allocate temporary solution matrix
    val Xtmp = DenseMatrix.zeros(math.max(A.numRows, A.numCols), nrhs);
    val M = if (!transpose) A.numRows else A.numCols;
    for (j <- 0 until nrhs; i <- 0 until M) { Xtmp(i,j) = V(i,j); }

    val newData = A.data.clone();

    // query optimal workspace
    val queryWork = new Array[Double](1);
    val queryInfo = new IntByReference(0);
    blas.dgels_(
      if (!transpose) "N" else "T",
      A.numRows, A.numCols, nrhs,
      newData, math.max(1,A.numRows),
      Xtmp.data, math.max(1,math.max(A.numRows,A.numCols)),
      queryWork, -1, queryInfo)
    
    // allocate workspace
    val work = {
      val lwork = {
        if (queryInfo.getValue != 0)
          math.max(1, math.min(A.numRows, A.numCols) + math.max(math.min(A.numRows, A.numCols), nrhs));
        else
          math.max(queryWork(0).toInt, 1);
      }
      new Array[Double](lwork);
    }

    // compute factorization
    val info = new IntByReference(0);
    blas.dgels_(
      if (!transpose) "N" else "T",
      A.numRows, A.numCols, nrhs,
      newData, math.max(1,A.numRows),
      Xtmp.data, math.max(1,math.max(A.numRows,A.numCols)),
      work, work.length, info);

    if (info.getValue< 0)
      throw new IllegalArgumentException;

    // extract solution
    val N = if (!transpose) A.numCols else A.numRows;
    for (j <- 0 until nrhs; i <- 0 until N) X(i,j) = Xtmp(i,j);

    X;
  }
}


class DenseMatrix(val numRows: Int, val numCols: Int, val data: Array[Double]) {
  require(numRows*numCols == data.size)
  
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
    val elems = Array.tabulate[Double](numCols)(j => this(i,j))
    new DenseMatrix(1, numCols, elems)
  }
  
  def apply(_slice: Slice, j: Int): DenseMatrix = {
    val elems = Array.tabulate[Double](numRows)(i => this(i,j))
    new DenseMatrix(numRows, 1, elems)
  }

  def update(i: Int, j: Int, x: Double) {
    data(index(i, j)) = x
  }
  
  def update(i: Int, _slice: Slice, x: DenseMatrix) {
    require(x.numRows == 1 && x.numCols == numCols)
    for (j <- 0 until numCols)
      this(i, j) = x(0, j)
  }

  def update(_slice: Slice, j: Int, x: DenseMatrix) {
    require(x.numRows == numRows && x.numCols == 1)
    for (i <- 0 until numRows)
      this(i, j) = x(i, 0)
  }
  
  def toComplex: DenseComplexMatrix = {
    val elems = new Array[Double](2*data.size)
    for (idx <- data.indices) {
      elems(2*idx+0) = data(idx)
      elems(2*idx+1) = 0
    }
    new DenseComplexMatrix(numRows, numCols, elems)
  }
  
  def transpose: DenseMatrix = {
    val ret = DenseMatrix.zeros(numCols, numRows)
    for (i <- 0 until numRows;
         j <- 0 until numCols) {
      ret(j, i) = this(i, j)
    }
    ret
  }
  
  def +(that: DenseMatrix): DenseMatrix = {
    require(numRows == that.numRows && numCols == that.numCols)
    val ret = DenseMatrix.zeros(numCols, numRows)
    for (i <- 0 until numRows;
         j <- 0 until numCols) {
      ret(i, j) = this(i, j) + that(i, j)
    }
    ret
  }

  def *(that: DenseMatrix): DenseMatrix = {
    require(numCols == that.numRows)
    val ret = DenseMatrix.zeros(numRows, that.numCols)
    val useBlas = true
    if (useBlas) {
      DenseMatrix.blasDgemm(ret, this, that)
    }
    else {
      for (i <- 0 until numRows;
           k <- 0 until numCols;
           j <- 0 until that.numCols) {
        ret(i, j) += this(i, k)*that(k, j)
      }
    }
    ret
  }
  
  def \(that: DenseMatrix): DenseMatrix = {
    require(numCols == that.numRows)
    require(that.numCols == 1)
    val ret = DenseMatrix.zeros(numRows, 1)
    DenseMatrix.QRSolve(ret, this, that, false)
    ret
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

