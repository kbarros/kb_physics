package kip.math.linalg4

import kip.math.Complex

// TODO: put implicits in object Dense (or somewhere lower priority)


trait DenseSlice
object :: extends DenseSlice


// Matrix types

trait Dense[S <: Scalar] extends Matrix[S, Dense] {
  val data: ScalarData[S]
  
  def index(i: Int, j: Int) = {
    checkKey(i, j)
    i + j*numRows // fortran column major convention
  }
  def indices = (for (j <- 0 until numCols; i <- 0 until numRows) yield (i, j)).toIterator
  
  def apply(i: Int, j: Int): S#A = data(index(i, j))
  def apply(i: Int, _slice: DenseSlice)(implicit mb: MatrixBuilder[S, DenseRow]): DenseRow[S] = {
    val ret = mb.zeros(1, numCols)
    for (j <- 0 until numCols) ret(0, j) = this(i, j)
    ret
  }
  def apply(_slice: DenseSlice, j: Int)(implicit mb: MatrixBuilder[S, DenseCol]): DenseCol[S] = {
    val ret = mb.zeros(numRows, 1)
    for (i <- 0 until numRows) ret(i, 0) = this(i, j)
    ret
  }
  
  def update(i: Int, j: Int, x: S#A) { data(index(i, j)) = x }
  def update[That[S <: Scalar] <: Matrix[S, That]](i: Int, _slice: DenseSlice, that: That[S]) {
    require(that.numRows == 1 && numCols == that.numCols, "Cannot perform matrix assignment: [%d, %d](%d, ::) <- [%d, %d]".format(
      numRows, numCols, i, that.numRows, that.numCols))
    for (j <- 0 until numCols) this(i, j) = that(0, j)
  }
  def update[That[S <: Scalar] <: Matrix[S, That]](_slice: DenseSlice, j: Int, that: That[S]) {
    require(that.numCols == 1 && numRows == that.numRows, "Cannot perform matrix assignment: [%d, %d](::, %d) <- [%d, %d]".format(
      numRows, numCols, j, that.numRows, that.numCols))
    for (i <- 0 until numRows) this(i, j) = that(i, 0)
  }

  def tranTo(that: Dense[S]): Unit = {
    require(numRows == that.numCols && numCols == that.numRows, "Cannot perform matrix transpose: [%d, %d] <- [%d, %d]^T".format(
      that.numRows, that.numCols, numRows, numCols))
    
    for (i <- 0 until math.min(numRows, numCols); j <- 0 to i) {
      val this_ij = this(i, j)
      val this_ji = this(j, i)
      that(i, j) = this_ji
      that(j, i) = this_ij
    }
    if (numCols > numRows) {
      for (i <- 0 until numRows; j <- numRows until numCols) {
        that(j, i) = this(i, j)
      }
    }
    if (numRows > numCols) {
      for (j <- 0 until numCols; i <- numCols until numRows) {
        that(j, i) = this(i, j)
      }
    }
  }
  
  def tran(implicit mb: MatrixBuilder[S, Dense]): Dense[S] = {
    val ret = mb.zeros(numCols, numRows)
    tranTo(ret)
    ret
  }
  
  def dag(implicit mb: MatrixBuilder[S, Dense]): Dense[S] = {
    val ret = tran
    ret.conjTo(ret)
    ret
  }
}

trait DenseRow[S <: Scalar] extends Matrix[S, DenseRow] with Dense[S] {
  def apply(i: Int): S#A = this(0, i)
  def update(i: Int, x: S#A): Unit = this(0, i) = x
  def tran(implicit mb: MatrixBuilder[S, DenseCol]): DenseCol[S] = {
    val ret = mb.zeros(numCols, numRows)
    tranTo(ret)
    ret
  }
  def dag(implicit mb: MatrixBuilder[S, DenseCol]): DenseCol[S] = {
    val ret = tran
    ret.conjTo(ret)
    ret
  }
  def *(that: DenseCol[S])
        (implicit mm: MatrixMultiplier[S, Dense, Dense, Dense],
                  mb: MatrixBuilder[S, Dense]): S#A = {
    val m = (this: Dense[S]) * (that: Dense[S])
    val ret = m(0, 0)
    m.data.dispose()
    ret
  }
}


trait DenseCol[S <: Scalar] extends Matrix[S, DenseCol] with Dense[S] {
  def apply(i: Int): S#A = this(i, 0)
  def update(i: Int, x: S#A): Unit = this(i, 0) = x
  def tran(implicit mb: MatrixBuilder[S, DenseRow]): DenseRow[S] = {
    val ret = mb.zeros(numCols, numRows)
    tranTo(ret)
    ret
  }
  def dag(implicit mb: MatrixBuilder[S, DenseRow]): DenseRow[S] = {
    val ret = tran
    ret.conjTo(ret)
    ret
  }
}

// Adders

trait DenseAdders {
  private def genericAddInPlace[S <: Scalar](sub: Boolean, m1: Dense[S], m2: Dense[S], ret: Dense[S]) {
    require(ret.numRows == m1.numRows && ret.numRows == m2.numRows,
            "Mismatched rows: %d, %d, %d".format(m1.numRows, m2.numRows, ret.numRows))
    require(ret.numCols == m1.numCols && ret.numCols == m2.numCols,
            "Mismatched cols: %d, %d, %d".format(m1.numCols, m2.numCols, ret.numCols))
    for ((i, j) <- ret.indices) ret(i, j) =
      if (sub) ret.scalar.sub(m1(i, j), m2(i, j)) else ret.scalar.add(m1(i, j), m2(i, j))
  }
  
  trait DDAdder[S <: Scalar] extends MatrixAdder[S, Dense, Dense, Dense] {
    def addInPlace(sub: Boolean, m1: Dense[S], m2: Dense[S], ret: Dense[S]) =
      genericAddInPlace(sub, m1, m2, ret)
  }
  implicit def dda[S <: Scalar] = new DDAdder[S] {}

  trait DRAdder[S <: Scalar] extends MatrixAdder[S, Dense, DenseRow, DenseRow] {
    def addInPlace(sub: Boolean, m1: Dense[S], m2: DenseRow[S], ret: DenseRow[S]) =
      genericAddInPlace(sub, m1, m2, ret)
  }
  implicit def dra[S <: Scalar] = new DRAdder[S] {}

  trait RDAdder[S <: Scalar] extends MatrixAdder[S, DenseRow, Dense, DenseRow] {
    def addInPlace(sub: Boolean, m1: DenseRow[S], m2: Dense[S], ret: DenseRow[S]) =
      genericAddInPlace(sub, m1, m2, ret)
  }
  implicit def rda[S <: Scalar] = new RDAdder[S] {}

  trait RRAdder[S <: Scalar] extends MatrixAdder[S, DenseRow, DenseRow, DenseRow] {
    def addInPlace(sub: Boolean, m1: DenseRow[S], m2: DenseRow[S], ret: DenseRow[S]) =
      genericAddInPlace(sub, m1, m2, ret)
  }
  implicit def rra[S <: Scalar] = new RRAdder[S] {}
}


// Multipliers

trait DenseMultipliers {
  private def genericGemm[S <: Scalar](alpha: S#A, beta: S#A, m1: Dense[S], m2: Dense[S], ret: Dense[S]) {
    require(ret.numRows == m1.numRows &&
            m1.numCols == m2.numRows &&
            m2.numCols == ret.numCols, "Cannot multiply matrices: (%d, %d) * (%d, %d) -> (%d, %d)".format(
              m1.numRows, m1.numCols, m2.numRows, m2.numCols, ret.numRows, ret.numCols))

//    if (Netlib.cblas == null) {
      for (i <- 0 until ret.numRows;
           k <- 0 until m1.numCols;
           j <- 0 until ret.numCols) {
        ret.data.madd(ret.index(i, j), m1.data, m1.index(i, k), m2.data, m2.index(k, j))
      }
//    }
/*
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
*/
  }
   
  trait DDMultiplier[S <: Scalar] extends MatrixMultiplier[S, Dense, Dense, Dense] {
    def gemm(alpha: S#A, beta: S#A, m1: Dense[S], m2: Dense[S], ret: Dense[S]) =
      genericGemm(alpha, beta, m1, m2, ret)
  }
  implicit def ddm[S <: Scalar] = new DDMultiplier[S] {}

  trait DCMultiplier[S <: Scalar] extends MatrixMultiplier[S, Dense, DenseCol, DenseCol] {
    def gemm(alpha: S#A, beta: S#A, m1: Dense[S], m2: DenseCol[S], ret: DenseCol[S]) =
      genericGemm(alpha, beta, m1, m2, ret)
  }
  implicit def dcm[S <: Scalar] = new DCMultiplier[S] {}

  trait RDMultiplier[S <: Scalar] extends MatrixMultiplier[S, DenseRow, Dense, DenseRow] {
    def gemm(alpha: S#A, beta: S#A, m1: DenseRow[S], m2: Dense[S], ret: DenseRow[S]) =
      genericGemm(alpha, beta, m1, m2, ret)
  }
  implicit def rdm[S <: Scalar] = new RDMultiplier[S] {}

  trait CRMultiplier[S <: Scalar] extends MatrixMultiplier[S, DenseCol, DenseRow, Dense] {
    def gemm(alpha: S#A, beta: S#A, m1: DenseCol[S], m2: DenseRow[S], ret: Dense[S]) =
      genericGemm(alpha, beta, m1, m2, ret)
  }
  implicit def crm[S <: Scalar] = new CRMultiplier[S] {}
}


// Builders
trait DenseBuilders {
  
  implicit def dense[S <: Scalar](implicit sb: ScalarData.Builder[S]) = new MatrixBuilder[S, Dense] {
    def zeros(numRows: Int, numCols: Int) = {
      require(numRows > 0 && numCols > 0, "Cannot build matrix with non-positive dimensions (%d, %d)".format(numRows, numCols))
      val nr = numRows
      val nc = numCols
      new Dense[S] {
        val data: ScalarData[S] = sb.build(nr*nc)
        val scalar: ScalarOps[S#A] = data.scalar
        val numRows = nr
        val numCols = nc
      }
    }
  }
  
  implicit def denseRow[S <: Scalar](implicit sb: ScalarData.Builder[S]) = new MatrixBuilder[S, DenseRow] {
    def zeros(numRows: Int, numCols: Int) = {
      require(numRows > 0 && numCols > 0, "Cannot build matrix with non-positive dimensions (%d, %d)".format(numRows, numCols))
      require(numRows == 1, "Cannot build row matrix with %d rows".format(numRows))
      val nc = numCols
      new DenseRow[S] {
        val data: ScalarData[S] = sb.build(1*nc)
        val scalar: ScalarOps[S#A] = data.scalar
        val numRows = 1
        val numCols = nc
      }
    }
  }
  
  implicit def denseCol[S <: Scalar](implicit sb: ScalarData.Builder[S]) = new MatrixBuilder[S, DenseCol] {
    def zeros(numRows: Int, numCols: Int) = {
      require(numRows > 0 && numCols > 0, "Cannot build matrix with non-positive dimensions (%d, %d)".format(numRows, numCols))
      require(numCols == 1, "Cannot build column matrix with %d cols".format(numCols))
      val nr = numRows
      new DenseCol[S] {
        val data: ScalarData[S] = sb.build(nr*1)
        val scalar: ScalarOps[S#A] = data.scalar
        val numRows = nr
        val numCols = 1
      }
    }
  }
  
  class DenseBuilderExtras[S <: Scalar](implicit sb: ScalarData.Builder[S]) {
    def zeros(numRows: Int, numCols: Int): Dense[S] = {
      dense.zeros(numRows, numCols)
    }

    def tabulate(numRows: Int, numCols: Int)(f: (Int, Int) => S#A): Dense[S] = {
      val m = dense.zeros(numRows, numCols)
      for ((i, j) <- m.indices) m(i, j) = f(i, j)
      m
    }

    def eye(numRows: Int): Dense[S] = {
      val m = dense.zeros(numRows, numRows)
      for (i <- 0 until numRows) m(i, i) = m.scalar.one
      m
    }
    
    def col(elems: S#A*): DenseCol[S] = {
      val m = denseCol.zeros(elems.size, 1)
      for (i <- 0 until m.numRows) m(i, 0) = elems(i)
      m
    }

    def row(elems: S#A*): DenseRow[S] = {
      val m = denseRow.zeros(1, elems.size)
      for (j <- 0 until m.numCols) m(0, j) = elems(j)
      m
    }

    def fromRows(row1: DenseRow[S], rows: DenseRow[S]*): Dense[S] = {
      require(row1.numRows == 1 && rows.forall(_.numRows == 1))
      require(rows.forall(_.numCols == row1.numCols))
      val ret = dense.zeros(1 + rows.size, row1.numCols)
      ret(0, ::) = row1
      for (i <- rows.indices)
        ret(i+1, ::) = rows(i)
      ret
    }
  }
  
  val denseRealFlt = new DenseBuilderExtras[Scalar.RealFlt]
  val denseRealDbl = new DenseBuilderExtras[Scalar.RealDbl]
  val denseComplexFlt = new DenseBuilderExtras[Scalar.ComplexFlt]
  val denseComplexDbl = new DenseBuilderExtras[Scalar.ComplexDbl]
}
