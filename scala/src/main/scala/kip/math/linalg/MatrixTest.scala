package kip.math
package linalg

/*
 * FUTURE WORK:
 *
 * - Sparse matrices
 * - CanBuildFrom for matrices
 * - JNA takes Buffer arguments
 * - ComplexDbl and ComplexFlt
 * 
 * */

object DenseMatrix {
  implicit object RealDbl    extends Builder[Double,  Array[Double]] ()(ScalarData.RealDbl, Netlib.RealDbl)
  implicit object RealFlt    extends Builder[Float,   Array[Float]]  ()(ScalarData.RealFlt, Netlib.RealFlt)
  implicit object ComplexDbl extends Builder[Complex, Array[Double]] ()(ScalarData.ComplexDbl, Netlib.ComplexDbl)
  implicit object ComplexFlt extends Builder[Complex, Array[Float]]  ()(ScalarData.ComplexFlt, Netlib.ComplexFlt)

  class Builder[@specialized(Float, Double) A, B](implicit scalar: ScalarData[A, B], netlib: Netlib[A, B]) {
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

  protected def extractReal(buffer: Any, i: Int): Double = {
    buffer match {
      case a: Array[Float] => a(i)
      case a: Array[Double] => a(i)
    }
  }

  /** (X := A \ V) or (X := A^T \ V) */
  def QRSolve[@specialized(Float, Double) A, B]
      (X: DenseMatrix[A, B], A: DenseMatrix[A, B], V: DenseMatrix[A, B], transpose: Boolean) = {
    require(X.numRows == A.numCols, "Wrong number of rows in return value");
    require(X.numCols == V.numCols, "Wrong number of rows in return value");
    val scalar = X.scalar
    val netlib = X.netlib
    val builder = X.matrix

    val nrhs = V.numCols;
    
    // allocate temporary solution matrix
    val Xtmp = builder.zeros(math.max(A.numRows, A.numCols), nrhs)
    val M = if (!transpose) A.numRows else A.numCols;
    for (j <- 0 until nrhs; i <- 0 until M) { Xtmp(i,j) = V(i,j) }

    val newData = scalar.clone(A.data)

    // query optimal workspace
    val queryWork = scalar.alloc(1);
    val queryInfo = new com.sun.jna.ptr.IntByReference(0);
    netlib.gels(
      if (!transpose) "N" else "T",
      A.numRows, A.numCols, nrhs,
      newData, math.max(1,A.numRows),
      Xtmp.data, math.max(1,math.max(A.numRows,A.numCols)),
      queryWork, -1, queryInfo)
    
    // allocate workspace
    val workSize =
      if (queryInfo.getValue != 0)
        math.max(1, math.min(A.numRows, A.numCols) + math.max(math.min(A.numRows, A.numCols), nrhs));
      else {
        math.max(extractReal(scalar.buffer(queryWork), 0).toInt, 1);
      }
    println("worksize "+workSize)
    val work = scalar.alloc(workSize);

    // compute factorization
    netlib.gels(
      if (!transpose) "N" else "T",
      A.numRows, A.numCols, nrhs,
      newData, math.max(1,A.numRows),
      Xtmp.data, math.max(1,math.max(A.numRows,A.numCols)),
      work, workSize, queryInfo);

    if (queryInfo.getValue< 0)
      throw new IllegalArgumentException;

    // extract solution
    val N = if (!transpose) A.numCols else A.numRows;
    for (j <- 0 until nrhs; i <- 0 until N) { X(i, j) = Xtmp(i, j) }

    X;
  }
}


class DenseMatrix[@specialized(Float, Double) A, B]
    (val numRows: Int, val numCols: Int, val data: B)
    (implicit val scalar: ScalarData[A, B], val netlib: Netlib[A, B]) {

  require(numRows*numCols == scalar.size(data))
  val matrix = new DenseMatrix.Builder
  
  override def clone(): DenseMatrix[A, B] = {
    new DenseMatrix(numRows, numCols, scalar.clone(data))
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
    if (Netlib.cblas == null) {
      for (i <- 0 until numRows;
           k <- 0 until numCols;
           j <- 0 until that.numCols) {
        scalar.madd(ret.data, ret.index(i, j), data, index(i, k), that.data, that.index(k, j))
      }
    }
    else {
      netlib.gemm(netlib.CblasColMajor, netlib.CblasNoTrans, netlib.CblasNoTrans,
                  ret.numRows, ret.numCols, // dimension of return matrix
                  numCols, // dimension of summation index
                  scalar.one, // alpha
                  data, numRows, // A matrix
                  that.data, that.numRows, // B matrix
                  scalar.zero, // beta
                  ret.data, ret.numRows // C matrix
                )
    }
    ret
  }
  
  def \(that: DenseMatrix[A, B]): DenseMatrix[A, B] = {
    require(numCols == that.numRows, "Size mismatch: cols %d != rows %d".format(numCols, numRows))
    require(that.numCols == 1)
    val ret = matrix.zeros(numRows, 1)
    DenseMatrix.QRSolve(ret, this, that, false)
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

