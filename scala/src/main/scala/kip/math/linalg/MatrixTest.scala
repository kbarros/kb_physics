package kip.math
package linalg

/*
 * FUTURE WORK:
 *
 * - LAPACK operations
 * - CanBuildFrom for matrices
 * - Sparse matrices
 * 
 * */

object DenseMatrix {
  implicit object RealDbl    extends Builder[Double,  Array[Double]] ()(ScalarData.RealDbl, NetlibOps.RealDbl)
  implicit object RealFlt    extends Builder[Float,   Array[Float]]  ()(ScalarData.RealFlt, NetlibOps.RealFlt)
  implicit object ComplexDbl extends Builder[Complex, Array[Double]] ()(ScalarData.ComplexDbl, NetlibOps.ComplexDbl)
  implicit object ComplexFlt extends Builder[Complex, Array[Float]]  ()(ScalarData.ComplexFlt, NetlibOps.ComplexFlt)

  class Builder[@specialized(Float, Double) A, B](implicit scalar: ScalarData[A, B], netlib: NetlibOps[A, B]) {
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
    (implicit scalar: ScalarData[A, B], netlib: NetlibOps[A, B]) {

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
    import CblasLib._
    netlib.gemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                ret.numRows, ret.numCols, // dimension of return matrix
                numCols, // dimension of summation index
                scalar.one, // alpha
                data, numRows, // A matrix
                that.data, that.numRows, // B matrix
                scalar.zero, // beta
                ret.data, ret.numRows // C matrix
              )
/*
    for (i <- 0 until numRows;
         k <- 0 until numCols;
         j <- 0 until that.numCols) {
      scalar.madd(ret.data, ret.index(i, j), data, index(i, k), that.data, that.index(k, j))
    }
*/
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


// ------------------

import kip.nativelib.VecLib
import com.sun.jna.ptr.IntByReference
import com.sun.jna.Native


object NetlibOps {
  implicit object RealFlt extends NetlibRealFltOps[Array[Float]]
  implicit object RealDbl extends NetlibRealDblOps[Array[Double]]
  implicit object ComplexFlt extends NetlibComplexFltOps[Array[Float]]
  implicit object ComplexDbl extends NetlibComplexDblOps[Array[Double]]
  
  lazy val cblas: CblasLib = Native.loadLibrary("vecLib", classOf[CblasLib]).asInstanceOf[CblasLib]
  lazy val lapack: LapackLib = Native.loadLibrary("vecLib", classOf[LapackLib]).asInstanceOf[LapackLib]
  
  object Implicits {
    implicit def intToIntByReference(a: Int) = new IntByReference(a)
    implicit def complexToDoubleArray(c: Complex) = Array[Double](c.re, c.im)
    implicit def complexToFloatArray(c: Complex) = Array[Float](c.re.toFloat, c.im.toFloat)
  }
}


trait NetlibOps[@specialized(Float, Double) A, B] {

  def gemm(Order: Int, TransA: Int, TransB: Int,
           M: Int, N: Int, K: Int,
           alpha: A,
           A: B, lda: Int,
           B: B, ldb: Int,
           beta: A,
           C: B, ldc: Int)

/*
  def gels[Buf](TRANS: String, M: Int, N: Int, NRHS: Int,
            A: Buf, LDA: Int,
            B: Buf, LDB: Int,
            WORK: Buf, LWORK: Int, INFO: IntByReference)

  def geev[Buf](JOBVL: String, JOBVR: String,
            N: Int, A: Buf, LDA: Int,
            WR: Buf, WI: Array[Double],
            VL: Buf, LDVL: Int,
            VR: Buf, LDVR: Int,
            WORK: Buf, LWORK: Int, INFO: IntByReference)
*/
}


class NetlibRealFltOps[B](implicit s: ScalarData[Float, B] with ScalarBufferedFlt) extends NetlibOps[Float, B] {
  import NetlibOps.cblas._
  
  def gemm(Order: Int, TransA: Int, TransB: Int,
           M: Int, N: Int, K: Int,
           alpha: Float,
           A: B, lda: Int,
           B: B, ldb: Int,
           beta: Float,
           C: B, ldc: Int) =
    cblas_sgemm(Order, TransA, TransB, M, N, K, alpha, s.buffer(A), lda, s.buffer(B), ldb, beta, s.buffer(C), ldc)
}


class NetlibRealDblOps[B](implicit s: ScalarData[Double, B] with ScalarBufferedDbl) extends NetlibOps[Double, B] {
  import NetlibOps.cblas._
  
  def gemm(Order: Int, TransA: Int, TransB: Int,
           M: Int, N: Int, K: Int,
           alpha: Double,
           A: B, lda: Int,
           B: B, ldb: Int,
           beta: Double,
           C: B, ldc: Int) =
    cblas_dgemm(Order, TransA, TransB, M, N, K, alpha, s.buffer(A), lda, s.buffer(B), ldb, beta, s.buffer(C), ldc)
}


class NetlibComplexFltOps[B](implicit s: ScalarData[Complex, B] with ScalarBufferedFlt) extends NetlibOps[Complex, B] {
  import NetlibOps.cblas._
  import NetlibOps.Implicits._
  
  def gemm(Order: Int, TransA: Int, TransB: Int,
           M: Int, N: Int, K: Int,
           alpha: Complex,
           A: B, lda: Int,
           B: B, ldb: Int,
           beta: Complex,
           C: B, ldc: Int) =
    cblas_cgemm(Order, TransA, TransB, M, N, K, alpha, s.buffer(A), lda, s.buffer(B), ldb, beta, s.buffer(C), ldc)
}


class NetlibComplexDblOps[B](implicit s: ScalarData[Complex, B] with ScalarBufferedDbl) extends NetlibOps[Complex, B] {
  import NetlibOps.cblas._
  import NetlibOps.Implicits._
  
  def gemm(Order: Int, TransA: Int, TransB: Int,
           M: Int, N: Int, K: Int,
           alpha: Complex,
           A: B, lda: Int,
           B: B, ldb: Int,
           beta: Complex,
           C: B, ldc: Int) =
    cblas_zgemm(Order, TransA, TransB, M, N, K, alpha, s.buffer(A), lda, s.buffer(B), ldb, beta, s.buffer(C), ldc)
}

