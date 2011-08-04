package kip.math.linalg4

import collection.mutable.HashMap


object Sparse extends SparseMultipliers with SparseBuilders

trait Sparse[S <: Scalar] extends Matrix[S, Sparse] {
  val data = new HashMap[(Int, Int), S#A]()
  
  override def apply(i: Int, j: Int): S#A = data.getOrElse((i, j), scalar.zero)
  override def update(i: Int, j: Int, x: S#A) { data((i, j)) = x }
  override def transform(f: S#A => S#A): this.type = {
    for ((i, j) <- data.keys) { this(i, j) = f(this(i, j)) }
    this
  }
  
  override def map[A2, S2 <: Scalar{type A=A2}, That[T <: Scalar] >: Sparse[T] <: Matrix[T, That]]
      (f: S#A => S2#A)(implicit scalar2: ScalarOps[S2], mb: MatrixBuilder[S2, That]): That[S2] = {
    require(f(scalar.zero) == scalar2.zero, "Map on sparse matrix must preserve zero")
    mb.map(this)(f)
  }

  def toDense(implicit mb: MatrixBuilder[S, Dense]): Dense[S] = {
    val ret = mb.zeros(numRows, numCols)
    for ((i, j) <- data.keys) { ret(i, j) = this(i, j) }
    ret
  }
}


trait SparseBuilders {
  implicit def sparseBuilder[S <: Scalar](implicit so: ScalarOps[S]) = new MatrixBuilder[S, Sparse] {
    def zeros(numRows: Int, numCols: Int) = {
      MatrixDims.checkDims(numRows, numCols)
      val nr = numRows
      val nc = numCols
      new Sparse[S] {
        override val scalar = so
        override val numRows = nr
        override val numCols = nc
      }
    }
    
    def duplicate(m: Sparse[S]): Sparse[S] = {
      val ret = zeros(m.numRows, m.numCols)
      for (k <- m.data.keys) { ret.data(k) = m.data(k) }
      ret
    }
    
    def transpose(m: Sparse[S]): Sparse[S] = {
      val ret = zeros(m.numCols, m.numRows)
      for ((i, j) <- m.data.keys) { ret.data((j, i)) = m.data((i, j)) }
      ret
    }
    
    def map[S0 <: Scalar](m: Sparse[S0])(f: S0#A => S#A): Sparse[S] = {
      val ret = zeros(m.numRows, m.numCols)
      for (k <- m.data.keys) { ret.data(k) = f(m.data(k)) }
      ret
    }
  }
}


trait SparseMultipliers {
  implicit def sparseDenseMultiplier[S <: Scalar] = new MatrixMultiplier[S, Sparse, Dense, Dense] {
    def gemm(alpha: S#A, beta: S#A, m1: Sparse[S], m2: Dense[S], ret: Dense[S]) {
      MatrixDims.checkMulTo(m1, m2, ret)
      ret.transform(_ => ret.scalar.zero)
      // ret_ij = \sum_k m1_ik m2_kj
      for ((i, k) <- m1.data.keys;
           val m1_ik = m1.apply(i, k); 
           j <- 0 until ret.numCols) {
        ret.scalar.madd(ret.data, ret.index(i, j), m2.data, m2.index(k, j), m1_ik)
      }
    }
  }
  implicit def denseSparseMultiplier[S <: Scalar] = new MatrixMultiplier[S, Dense, Sparse, Dense] {
    def gemm(alpha: S#A, beta: S#A, m1: Dense[S], m2: Sparse[S], ret: Dense[S]) {
      MatrixDims.checkMulTo(m1, m2, ret)
      ret.transform(_ => ret.scalar.zero)
      // ret_ij = \sum_k m1_ik m2_kj
      for ((k, j) <- m2.data.keys;
           val m2_kj = m1.apply(k, j); 
           i <- 0 until ret.numRows) {
        ret.scalar.madd(ret.data, ret.index(i, j), m1.data, m1.index(i, k), m2_kj)
      }
    }
  }
}
