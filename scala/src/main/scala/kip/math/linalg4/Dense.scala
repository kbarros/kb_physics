package kip.math.linalg4

import kip.math.Complex


trait DenseSlice
object :: extends DenseSlice


object Dense {
  implicit def toDenseRealLapackOps   [S <: Scalar.RealTyp]   (m: Dense[S]) = new DenseRealLapackOps(m)
  implicit def toDenseComplexLapackOps[S <: Scalar.ComplexTyp](m: Dense[S]) = new DenseComplexLapackOps(m)
}


trait Dense[S <: Scalar] extends Matrix[S, Dense] {
  val netlib: Netlib[S]
  val data: RawData[S#Raw, S#Buf]
  
  def index(i: Int, j: Int) = {
    checkKey(i, j)
    i + j*numRows // fortran column major convention
  }
  def indices = for (j <- (0 until numCols).view; i <- (0 until numRows).view) yield (i, j)
  
  def apply(i: Int, j: Int): S#A = scalar.read(data, index(i, j))
  def apply(i: Int): S#A = {
    if (numRows == 1)
      this(0, i)
    else if (numCols == 1)
      this(i, 0)
    else {
      require(numRows == 1 || numCols == 1, "Cannot apply a single index to non-vector matrix of shape [%d, %d]".format(numRows, numCols))
      sys.error("")
    }
  }
  def apply(i: Int, _slice: DenseSlice)(implicit mb: MatrixBuilder[S, Dense]): Dense[S] = {
    val ret = mb.zeros(1, numCols)
    for (j <- 0 until numCols) ret(0, j) = this(i, j)
    ret
  }
  def apply(_slice: DenseSlice, j: Int)(implicit mb: MatrixBuilder[S, Dense]): Dense[S] = {
    val ret = mb.zeros(numRows, 1)
    for (i <- 0 until numRows) ret(i, 0) = this(i, j)
    ret
  }
  
  def update(i: Int, j: Int, x: S#A): Unit = scalar.write(data, index(i, j), x)
  def update(i: Int, x: S#A): Unit = {
    if (numRows == 1)
      this(0, i) = x
    else if (numCols == 1)
      this(i, 0) = x
    else {
      require(numRows == 1 || numCols == 1, "Cannot apply a single index to non-vector matrix of shape [%d, %d]".format(numRows, numCols))
      sys.error("")
    }
  }
  def update[That[S <: Scalar] <: Matrix[S, That]](i: Int, _slice: DenseSlice, that: That[S]) {
    require(that.numRows == 1 && numCols == that.numCols, "Cannot perform matrix assignment of shape: [%d, %d] -> [%d, %d](%d, ::)".format(
      that.numRows, that.numCols, numRows, numCols, i))
    for (j <- 0 until numCols) this(i, j) = that(0, j)
  }
  def update[That[S <: Scalar] <: Matrix[S, That]](_slice: DenseSlice, j: Int, that: That[S]) {
    require(that.numCols == 1 && numRows == that.numRows, "Cannot perform matrix assignment of shape: [%d, %d] -> [%d, %d](::, %d)".format(
      that.numRows, that.numCols, numRows, numCols, j))
    for (i <- 0 until numRows) this(i, j) = that(i, 0)
  }
  
  def tran(implicit mb: MatrixBuilder[S, Dense]): Dense[S] = {
    val ret = mb.zeros(numCols, numRows)
    for (i <- 0 until math.min(numRows, numCols); j <- 0 to i) {
      val this_ij = this(i, j)
      val this_ji = this(j, i)
      ret(i, j) = this_ji
      ret(j, i) = this_ij
    }
    if (numCols > numRows) {
      for (i <- 0 until numRows; j <- numRows until numCols) {
        ret(j, i) = this(i, j)
      }
    }
    if (numRows > numCols) {
      for (j <- 0 until numCols; i <- numCols until numRows) {
        ret(j, i) = this(i, j)
      }
    }
    ret
  }
  
  def dag(implicit mb: MatrixBuilder[S, Dense]): Dense[S] = {
    tran.transform(scalar.conj(_))
  }

  def dot(that: Dense[S])(implicit mm: MatrixMultiplier[S, Dense, Dense, Dense],
      mb: MatrixBuilder[S, Dense]): S#A = {
    val m = (this: Dense[S]) * (that: Dense[S])
    val ret = m(0, 0)
    m.data.dispose()
    ret
  }
}


trait DenseAdders {
  trait DenseDenseAdder[S <: Scalar] extends MatrixAdder[S, Dense, Dense, Dense] {
    def addInPlace(sub: Boolean, m1: Dense[S], m2: Dense[S], ret: Dense[S]) = {
      require(
          ret.numRows == m1.numRows &&
          ret.numRows == m2.numRows &&
          ret.numCols == m1.numCols &&
          ret.numCols == m2.numCols, "Cannot add matrices of shape: [%d, %d] + [%d, %d] -> [%d, %d]".format(
              m1.numRows, m1.numCols, m2.numRows, m2.numCols, ret.numRows, ret.numCols
          )
      )
      for ((i, j) <- ret.indices) ret(i, j) =
        if (sub) ret.scalar.sub(m1(i, j), m2(i, j)) else ret.scalar.add(m1(i, j), m2(i, j))
    }
  }
  implicit def denseDenseAdder[S <: Scalar] = new DenseDenseAdder[S] {}
}


trait DenseMultipliers {
  trait DenseDenseMultiplier[S <: Scalar] extends MatrixMultiplier[S, Dense, Dense, Dense] {
    def gemm(alpha: S#A, beta: S#A, m1: Dense[S], m2: Dense[S], ret: Dense[S]) {
      require(
          ret.numRows == m1.numRows &&
          m1.numCols == m2.numRows &&
          m2.numCols == ret.numCols, "Cannot multiply matrices of shape: [%d, %d] * [%d, %d] -> [%d, %d]".format(
              m1.numRows, m1.numCols, m2.numRows, m2.numCols, ret.numRows, ret.numCols
          )
      )

      if (ret.netlib == null) {
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
  }
  implicit def denseDenseMultiplier[S <: Scalar] = new DenseDenseMultiplier[S] {}
}



trait DenseBuilders {
  
  class DenseBuilder[S <: Scalar](implicit so: ScalarOps[S], sb: RawData.Builder[S#Raw, S#Buf], nl: Netlib[S]) extends MatrixBuilder[S, Dense] {
    def zeros(numRows: Int, numCols: Int) = {
      require(numRows > 0 && numCols > 0, "Cannot build matrix with non-positive dimensions [%d, %d]".format(numRows, numCols))
      val nr = numRows
      val nc = numCols
      new Dense[S] {
        val netlib: Netlib[S] = nl
        val data: RawData[S#Raw, S#Buf] = sb.build(so.components*nr*nc)
        // TODO remove explicit types
        val scalar: ScalarOps[S] = so
        val numRows = nr
        val numCols = nc
      }
    }

    def tabulate(numRows: Int, numCols: Int)(f: (Int, Int) => S#A): Dense[S] = {
      val m = zeros(numRows, numCols)
      for ((i, j) <- m.indices) m(i, j) = f(i, j)
      m
    }

    def eye(numRows: Int): Dense[S] = {
      val m = zeros(numRows, numRows)
      for (i <- 0 until numRows) m(i, i) = m.scalar.one
      m
    }
    
    def col(elems: S#A*): Dense[S] = {
      val m = zeros(elems.size, 1)
      for (i <- 0 until m.numRows) m(i, 0) = elems(i)
      m
    }

    def row(elems: S#A*): Dense[S] = {
      val m = zeros(1, elems.size)
      for (j <- 0 until m.numCols) m(0, j) = elems(j)
      m
    }

    def fromRows(row1: Dense[S], rows: Dense[S]*): Dense[S] = {
      require(row1.numRows == 1 && rows.forall(_.numRows == 1))
      require(rows.forall(_.numCols == row1.numCols))
      val ret = zeros(1 + rows.size, row1.numCols)
      ret(0, ::) = row1
      for (i <- rows.indices)
        ret(i+1, ::) = rows(i)
      ret
    }
  }
  
  implicit def denseBuilder[S <: Scalar](implicit so: ScalarOps[S], sb: RawData.Builder[S#Raw, S#Buf], nl: Netlib[S]) = new DenseBuilder
  
  val denseRealFlt    = new DenseBuilder[Scalar.RealFlt]
  val denseRealDbl    = new DenseBuilder[Scalar.RealDbl]
  val denseComplexFlt = new DenseBuilder[Scalar.ComplexFlt]
  val denseComplexDbl = new DenseBuilder[Scalar.ComplexDbl]
}


// --------------------------------------
// Lapack operations


object DenseLapackOps {
    /** A \ V -> X where A X ~= V */
  def QRSolve[S <: Scalar](X: Dense[S], A: Dense[S], V: Dense[S], transpose: Boolean)(implicit mb: MatrixBuilder[S, Dense]) {
    require(
        A.numCols == X.numRows && A.numRows == V.numRows && X.numCols == V.numCols, 
        "Cannot divide matrices: [%d, %d] \\ [%d, %d] -> [%d, %d]".format(
            A.numRows, A.numCols, V.numRows, V.numCols, X.numRows, X.numCols 
        )
    )
    require(X.netlib != null, "Netlib library required for division operation")
    
    val nrhs = V.numCols;
    
    // allocate temporary solution matrix
    // TODO: is X big enough to use directly?
    val Xtmp = mb.zeros(math.max(A.numRows, A.numCols), nrhs)
    val M = if (!transpose) A.numRows else A.numCols;
    for (j <- 0 until nrhs; i <- 0 until M) { Xtmp(i, j) = V(i, j) }

    val newData = A.duplicate(mb);

    // query optimal workspace
    val queryWork = mb.zeros(1, 1)
    val queryInfo = new com.sun.jna.ptr.IntByReference(0);
    X.netlib.gels(
      if (!transpose) "N" else "T",
      A.numRows, A.numCols, nrhs,
      newData.data.buffer, math.max(1,A.numRows),
      Xtmp.data.buffer, math.max(1,math.max(A.numRows,A.numCols)),
      queryWork.data.buffer, -1, queryInfo)
    
    // allocate workspace
    val workSize =
      if (queryInfo.getValue != 0)
        math.max(1, math.min(A.numRows, A.numCols) + math.max(math.min(A.numRows, A.numCols), nrhs));
      else
        math.max(1, X.netlib.readBufferInt(queryWork.data.buffer, 0));
    val work = mb.zeros(workSize, 1)

    // compute factorization
    X.netlib.gels(
      if (!transpose) "N" else "T",
      A.numRows, A.numCols, nrhs,
      newData.data.buffer, math.max(1,A.numRows),
      Xtmp.data.buffer, math.max(1,math.max(A.numRows,A.numCols)),
      work.data.buffer, workSize, queryInfo);

    if (queryInfo.getValue< 0)
      throw new IllegalArgumentException;

    // extract solution
    val N = if (!transpose) A.numCols else A.numRows;
    for (j <- 0 until nrhs; i <- 0 until N) X(i, j) = Xtmp(i, j)
  }
}


class DenseLapackOps[S <: Scalar](self: Dense[S]) {
  def \ (V: Dense[S])(implicit mb: MatrixBuilder[S, Dense]): Dense[S] = {
    val ret = mb.zeros(V.numRows, V.numCols)
    DenseLapackOps.QRSolve(ret, self, V, false)
    ret
  }
}

class DenseRealLapackOps[S <: Scalar.RealTyp](self: Dense[S]) extends DenseLapackOps(self) {
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


class DenseComplexLapackOps[S <: Scalar.ComplexTyp](self: Dense[S]) extends DenseLapackOps(self) {
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
