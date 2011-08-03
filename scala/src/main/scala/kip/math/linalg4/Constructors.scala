package kip.math.linalg4


object Constructors {
  val realFlt    = new Constructors[Scalar.RealFlt]
  val realDbl    = new Constructors[Scalar.RealDbl]
  val complexFlt = new Constructors[Scalar.ComplexFlt]
  val complexDbl = new Constructors[Scalar.ComplexDbl]
}


class Constructors[S <: Scalar]() {
  // --------------------------------------
  // Dense constructors
  
  def dense(numRows: Int, numCols: Int)(implicit mb: MatrixBuilder[S, Dense]): Dense[S] = {
    mb.zeros(numRows, numCols)
  }

  def tabulate(numRows: Int, numCols: Int)(f: (Int, Int) => S#A)(implicit mb: MatrixBuilder[S, Dense]): Dense[S] = {
      val ret = dense(numRows, numCols)
      for (i <- 0 until numRows; j <- 0 until numCols) { ret(i, j) = f(i, j) }
      ret
  }

  def eye(numRows: Int)(implicit mb: MatrixBuilder[S, Dense]): Dense[S] = {
      val m = dense(numRows, numRows)
      for (i <- 0 until numRows) m(i, i) = m.scalar.one
      m
  }

  def col(elems: S#A*)(implicit mb: MatrixBuilder[S, Dense]): Dense[S] = {
      val m = dense(elems.size, 1)
      for (i <- 0 until m.numRows) m(i, 0) = elems(i)
      m
  }

  def row(elems: S#A*)(implicit mb: MatrixBuilder[S, Dense]): Dense[S] = {
      val m = dense(1, elems.size)
      for (j <- 0 until m.numCols) m(0, j) = elems(j)
      m
  }

  def fromRows(row1: Dense[S], rows: Dense[S]*)(implicit mb: MatrixBuilder[S, Dense]): Dense[S] = {
      require(row1.numRows == 1 && rows.forall(_.numRows == 1))
      require(rows.forall(_.numCols == row1.numCols))
      val ret = dense(1 + rows.size, row1.numCols)
      ret(0, ::) = row1
      for (i <- rows.indices)
        ret(i+1, ::) = rows(i)
        ret
  }
  
  // --------------------------------------
  // Sparse constructors
  
  def sparse(numRows: Int, numCols: Int)(implicit mb: MatrixBuilder[S, Sparse]): Sparse[S] = {
    mb.zeros(numRows, numCols)
  }

}

