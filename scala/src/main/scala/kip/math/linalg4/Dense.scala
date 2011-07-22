package kip.math.linalg4

import kip.math.Complex


trait DenseSlice
object :: extends DenseSlice


object Dense {
  implicit def toDenseRealOps   [S <: Scalar.RealTyp]   (m: Dense[S]) = new DenseRealOps(m)
  implicit def toDenseComplexOps[S <: Scalar.ComplexTyp](m: Dense[S]) = new DenseComplexOps(m)
}


// Matrix types

trait Dense[S <: Scalar] extends Matrix[S, Dense] {
  val netlib: Netlib[S]
  
  val data: RawData[S#Raw, S#Buf]
  
  def index(i: Int, j: Int) = {
    checkKey(i, j)
    i + j*numRows // fortran column major convention
  }
  def indices = (for (j <- 0 until numCols; i <- 0 until numRows) yield (i, j)).toIterator
  
  def apply(i: Int, j: Int): S#A = scalar.read(data, index(i, j))
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
  
  def update(i: Int, j: Int, x: S#A): Unit = scalar.write(data, index(i, j), x)
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

  trait CCAdder[S <: Scalar] extends MatrixAdder[S, DenseCol, DenseCol, DenseCol] {
    def addInPlace(sub: Boolean, m1: DenseCol[S], m2: DenseCol[S], ret: DenseCol[S]) =
      genericAddInPlace(sub, m1, m2, ret)
  }
  implicit def cca[S <: Scalar] = new CCAdder[S] {}
}


// Multipliers

trait DenseMultipliers {
  private def genericGemm[S <: Scalar](alpha: S#A, beta: S#A, m1: Dense[S], m2: Dense[S], ret: Dense[S]) {
    require(ret.numRows == m1.numRows &&
            m1.numCols == m2.numRows &&
            m2.numCols == ret.numCols, "Cannot multiply matrices: [%d, %d] * [%d, %d] -> [%d, %d]".format(
              m1.numRows, m1.numCols, m2.numRows, m2.numCols, ret.numRows, ret.numCols))

    if (Netlib.cblas == null) {
      for (i <- 0 until ret.numRows;
           k <- 0 until m1.numCols;
           j <- 0 until ret.numCols) {
        ret.scalar.madd(ret.data, ret.index(i, j), m1.data, m1.index(i, k), m2.data, m2.index(k, j))
      }
    }
    else {
      ret.netlib.gemm(Netlib.CblasColMajor, Netlib.CblasNoTrans, Netlib.CblasNoTrans,
                      ret.numRows, ret.numCols, // dimension of return matrix
                      m1.numCols, // dimension of summation index
                      ret.scalar.one, // alpha
                      m1.data.buffer, m1.numRows, // A matrix
                      m2.data.buffer, m2.numRows, // B matrix
                      ret.scalar.zero, // beta
                      ret.data.buffer, ret.numRows // C matrix
                    )
    }
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
  
  implicit def dense[S <: Scalar](implicit so: ScalarOps[S], sb: RawData.Builder[S#Raw, S#Buf], nl: Netlib[S]) = new MatrixBuilder[S, Dense] {
    def zeros(numRows: Int, numCols: Int) = {
      require(numRows > 0 && numCols > 0, "Cannot build matrix with non-positive dimensions [%d, %d]".format(numRows, numCols))
      val nr = numRows
      val nc = numCols
      new Dense[S] {
        val netlib: Netlib[S] = nl
        val data: RawData[S#Raw, S#Buf] = sb.build(nr*nc)
        // TODO remove explicit types
        val scalar: ScalarOps[S] = so
        val numRows = nr
        val numCols = nc
      }
    }
  }
  
  // TODO: remove row, col
  implicit def denseRow[S <: Scalar](implicit so: ScalarOps[S], sb: RawData.Builder[S#Raw, S#Buf], nl: Netlib[S]) = new MatrixBuilder[S, DenseRow] {
    def zeros(numRows: Int, numCols: Int) = {
      require(numRows > 0 && numCols > 0, "Cannot build matrix with non-positive dimensions [%d, %d]".format(numRows, numCols))
      require(numRows == 1, "Cannot build row matrix with %d rows".format(numRows))
      val nc = numCols
      new DenseRow[S] {
        val netlib: Netlib[S] = nl
        val data: RawData[S#Raw, S#Buf] = sb.build(1*nc)
        val scalar: ScalarOps[S] = so
        val numRows = 1
        val numCols = nc
      }
    }
  }
  
  implicit def denseCol[S <: Scalar](implicit so: ScalarOps[S], sb: RawData.Builder[S#Raw, S#Buf], nl: Netlib[S]) = new MatrixBuilder[S, DenseCol] {
    def zeros(numRows: Int, numCols: Int) = {
      require(numRows > 0 && numCols > 0, "Cannot build matrix with non-positive dimensions [%d, %d]".format(numRows, numCols))
      require(numCols == 1, "Cannot build column matrix with %d cols".format(numCols))
      val nr = numRows
      new DenseCol[S] {
        val netlib: Netlib[S] = nl
        val data: RawData[S#Raw, S#Buf] = sb.build(nr*1)
        val scalar: ScalarOps[S] = so
        val numRows = nr
        val numCols = 1
      }
    }
  }
  
  class DenseBuilderExtras[S <: Scalar](implicit so: ScalarOps[S], sb: RawData.Builder[S#Raw, S#Buf], nl: Netlib[S]) {
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
  
  val denseRealFlt    = new DenseBuilderExtras[Scalar.RealFlt]
  val denseRealDbl    = new DenseBuilderExtras[Scalar.RealDbl]
  val denseComplexFlt = new DenseBuilderExtras[Scalar.ComplexFlt]
  val denseComplexDbl = new DenseBuilderExtras[Scalar.ComplexDbl]
}

// Lapack operations


class DenseRealOps[S <: Scalar.RealTyp](self: Dense[S]) {
  def eig(implicit mb: MatrixBuilder[S, Dense]): (Dense[S], Dense[S], Dense[S]) = {
    self.netlib match {
      case netlib: NetlibReal[_] => {
        require(self.numRows == self.numCols, "Cannot find eigenvectors of non-square matrix [%d, %d]".format(self.numRows, self.numCols))
        
        val n = self.numRows

        // Allocate space for the decomposition
        val Wr = mb.zeros(n, 1)
        val Wi = mb.zeros(n, 1)
        val Vr = mb.zeros(n, n)
        
        // Find the needed workspace
        val worksize = mb.zeros(1, 1)
        val info = new com.sun.jna.ptr.IntByReference(0)
        netlib.geev(
          "N", "V", // compute right eigenvectors only
          n,
          netlib.emptyBuffer, math.max(1,n),
          netlib.emptyBuffer, netlib.emptyBuffer,
          netlib.emptyBuffer, math.max(1,n),
          netlib.emptyBuffer, math.max(1,n),
          worksize.data.buffer, -1, info)

        // Allocate the workspace
        val lwork: Int =
          if (info.getValue != 0)
            math.max(1, 4*n)
          else
            math.max(1, netlib.readBufferInt(worksize.data.buffer, 0))
        val work = mb.zeros(lwork, 1)

        // Factor it!
        netlib.geev(
          "N", "V", n,
          self.duplicate.data.buffer, math.max(1,n),
          Wr.data.buffer, Wi.data.buffer,
          netlib.emptyBuffer, math.max(1, n),
          Vr.data.buffer, math.max(1,n),
          work.data.buffer, lwork, info)

        require(info.getValue >= 0, "Error in dgeev argument %d".format(-info.getValue))
        require(info.getValue <= 0, "Not converged dgeev; only %d of %d eigenvalues computed".format(info.getValue, self.numRows))
        
        (Wr, Wi, Vr)
      }
      case _ => sys.error("Netlib sgeev/dgeev operation unavailable.")
    }
  }
}


class DenseComplexOps[S <: Scalar.ComplexTyp](self: Dense[S]) {
  def eig(implicit mb: MatrixBuilder[S, Dense]): (Dense[S], Dense[S]) = {
    self.netlib match {
      case netlib: NetlibComplex[_] => {
        require(self.numRows == self.numCols, "Cannot find eigenvectors of non-square matrix [%d, %d]".format(self.numRows, self.numCols))
        
        val n = self.numRows

        // Allocate space for the decomposition
        val W = mb.zeros(n, 1)
        val Vr = mb.zeros(n, n)

        // Find the needed workspace
        val worksize = mb.zeros(1, 1)
        val info = new com.sun.jna.ptr.IntByReference(0)
        netlib.geev(
          "N", "V", // compute right eigenvectors only
          n,
          netlib.emptyBuffer, math.max(1,n),
          netlib.emptyBuffer,
          netlib.emptyBuffer, math.max(1,n),
          netlib.emptyBuffer, math.max(1,n),
          worksize.data.buffer, -1, netlib.emptyBuffer, info)

        // Allocate the workspace
        val lwork: Int =
          if (info.getValue != 0)
            math.max(1, 4*n)
          else
            math.max(1, netlib.readBufferInt(worksize.data.buffer, 0))
        val work = mb.zeros(lwork, 1)
        val rwork = mb.zeros(n, 1) // 2*n float/double values

        // Factor it!
        netlib.geev(
          "N", "V", n,
          self.duplicate.data.buffer, math.max(1,n),
          W.data.buffer,
          netlib.emptyBuffer, math.max(1, n),
          Vr.data.buffer, math.max(1,n),
          work.data.buffer, lwork, rwork.data.buffer, info)

        require(info.getValue >= 0, "Error in dgeev argument %d".format(-info.getValue))
        require(info.getValue <= 0, "Not converged dgeev; only %d of %d eigenvalues computed".format(info.getValue, self.numRows))
        
        (W, Vr)
      }
      case _ => sys.error("Netlib cgeev/zgeev operation unavailable.")
    }
  }
}
