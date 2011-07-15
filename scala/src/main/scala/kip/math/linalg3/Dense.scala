package kip.math.linalg3

import kip.math.Complex
import kip.math.linalg2.{ScalarOps, ScalarData}


// TODO: put implicits in object Dense (or somewhere lower priority)


trait DenseSlice
object :: extends DenseSlice


// Matrix types

trait Dense[A, Raw] extends Matrix[A, Raw, Dense] {
  val data: ScalarData[A, Raw]
  
  def index(i: Int, j: Int) = {
    checkKey(i, j)
    i + j*numRows // fortran column major convention
  }
  def indices = (for (j <- 0 until numCols; i <- 0 until numRows) yield (i, j)).toIterator
  
  
  def apply(i: Int, j: Int): A = data(index(i, j))
  def apply(i: Int, _slice: DenseSlice)(implicit mb: MatrixBuilder[A, Raw, DenseRow]): DenseRow[A, Raw] = {
    val ret = mb.zeros(1, numCols)
    for (j <- 0 until numCols) ret(0, j) = this(i, j)
    ret
  }
  def apply(_slice: DenseSlice, j: Int)(implicit mb: MatrixBuilder[A, Raw, DenseCol]): DenseCol[A, Raw] = {
    val ret = mb.zeros(numRows, 1)
    for (i <- 0 until numRows) ret(i, 0) = this(i, j)
    ret
  }

  def update(i: Int, j: Int, x: A) { data(index(i, j)) = x }
  def update[That[C,D] <: Matrix[C, D, That]](i: Int, _slice: DenseSlice, that: That[A, Raw]) {
    require(that.numRows == 1 && numCols == that.numCols, "Cannot perform matrix assignment: [%d, %d](%d, ::) <- [%d, %d]".format(
      numRows, numCols, i, that.numRows, that.numCols))
    for (j <- 0 until numCols) this(i, j) = that(0, j)
  }
  def update[That[C,D] <: Matrix[C, D, That]](_slice: DenseSlice, j: Int, that: That[A, Raw]) {
    require(that.numCols == 1 && numRows == that.numRows, "Cannot perform matrix assignment: [%d, %d](::, %d) <- [%d, %d]".format(
      numRows, numCols, j, that.numRows, that.numCols))
    for (i <- 0 until numRows) this(i, j) = that(i, 0)
  }
  
  def tran(implicit mb: MatrixBuilder[A, Raw, Dense]): Dense[A, Raw] = {
    val ret = mb.zeros(numCols, numRows)
    for (i <- 0 until ret.numRows; j <- 0 until ret.numCols) ret(i, j) = this(j, i)
    ret
  }
  
  def dag(implicit mb: MatrixBuilder[A, Raw, Dense]): Dense[A, Raw] = {
    val ret = tran
    ret.conjTo(ret)
    ret
  }

  def tabulate(f: (Int, Int) => A) {
    for ((i, j) <- indices) this(i, j) = f(i, j)
  }
}

trait DenseRow[A, Raw] extends Matrix[A, Raw, DenseRow] with Dense[A, Raw] {
  def apply(i: Int): A = this(0, i)
  def update(i: Int, x: A): Unit = this(0, i) = x
  def *(that: DenseCol[A, Raw]): A = {
    scalar.zero
  }
  def tran(implicit mb: MatrixBuilder[A, Raw, DenseCol]): DenseCol[A, Raw] = {
    val ret = mb.zeros(numCols, numRows)
    for (i <- 0 until ret.numRows; j <- 0 until ret.numCols) ret(i, j) = this(j, i)
    ret
  }
  def dag(implicit mb: MatrixBuilder[A, Raw, DenseCol]): DenseCol[A, Raw] = {
    val ret = tran
    ret.conjTo(ret)
    ret
  }
}

trait DenseCol[A, Raw] extends Matrix[A, Raw, DenseCol] with Dense[A, Raw] {
  def apply(i: Int): A = this(i, 0)
  def update(i: Int, x: A): Unit = this(i, 0) = x
  def tran(implicit mb: MatrixBuilder[A, Raw, DenseRow]): DenseRow[A, Raw] = {
    val ret = mb.zeros(numCols, numRows)
    for (i <- 0 until ret.numRows; j <- 0 until ret.numCols) ret(i, j) = this(j, i)
    ret
  }
  def dag(implicit mb: MatrixBuilder[A, Raw, DenseRow]): DenseRow[A, Raw] = {
    val ret = tran
    ret.conjTo(ret)
    ret
  }
}


// Adders

trait DenseAdders {
  private def genericAddInPlace[A, Raw](sub: Boolean, m1: Dense[A, Raw], m2: Dense[A, Raw], ret: Dense[A, Raw]) {
    require(ret.numRows == m1.numRows && ret.numRows == m2.numRows,
            "Mismatched rows: %d, %d, %d".format(m1.numRows, m2.numRows, ret.numRows))
    require(ret.numCols == m1.numCols && ret.numCols == m2.numCols,
            "Mismatched cols: %d, %d, %d".format(m1.numCols, m2.numCols, ret.numCols))
    for ((i, j) <- ret.indices) ret(i, j) =
      if (sub) ret.scalar.sub(m1(i, j), m2(i, j)) else ret.scalar.add(m1(i, j), m2(i, j))
  }
  
  trait DDAdder[A, Raw] extends MatrixAdder[A, Raw, Dense, Dense, Dense] {
    def addInPlace(sub: Boolean, m1: Dense[A, Raw], m2: Dense[A, Raw], ret: Dense[A, Raw]) =
      genericAddInPlace(sub, m1, m2, ret)
  }
  implicit def dda[A, Raw] = new DDAdder[A, Raw] {}

  trait DRAdder[A, Raw] extends MatrixAdder[A, Raw, Dense, DenseRow, DenseRow] {
    def addInPlace(sub: Boolean, m1: Dense[A, Raw], m2: DenseRow[A, Raw], ret: DenseRow[A, Raw]) =
      genericAddInPlace(sub, m1, m2, ret)
  }
  implicit def dra[A, Raw] = new DRAdder[A, Raw] {}

  trait RDAdder[A, Raw] extends MatrixAdder[A, Raw, DenseRow, Dense, DenseRow] {
    def addInPlace(sub: Boolean, m1: DenseRow[A, Raw], m2: Dense[A, Raw], ret: DenseRow[A, Raw]) =
      genericAddInPlace(sub, m1, m2, ret)
  }
  implicit def rda[A, Raw] = new RDAdder[A, Raw] {}

  trait RRAdder[A, Raw] extends MatrixAdder[A, Raw, DenseRow, DenseRow, DenseRow] {
    def addInPlace(sub: Boolean, m1: DenseRow[A, Raw], m2: DenseRow[A, Raw], ret: DenseRow[A, Raw]) =
      genericAddInPlace(sub, m1, m2, ret)
  }
  implicit def rra[A, Raw] = new RRAdder[A, Raw] {}
}


// Multipliers

trait DenseMultipliers {
  private def genericGemm[A, Raw](alpha: A, beta: A, m1: Dense[A, Raw], m2: Dense[A, Raw], ret: Dense[A, Raw]) {
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
   
  trait DDMultiplier[A, Raw] extends MatrixMultiplier[A, Raw, Dense, Dense, Dense] {
    def gemm(alpha: A, beta: A, m1: Dense[A, Raw], m2: Dense[A, Raw], ret: Dense[A, Raw]) =
      genericGemm(alpha, beta, m1, m2, ret)
  }
  implicit def ddm[A, Raw] = new DDMultiplier[A, Raw] {}

  trait DCMultiplier[A, Raw] extends MatrixMultiplier[A, Raw, Dense, DenseCol, DenseCol] {
    def gemm(alpha: A, beta: A, m1: Dense[A, Raw], m2: DenseCol[A, Raw], ret: DenseCol[A, Raw]) =
      genericGemm(alpha, beta, m1, m2, ret)
  }
  implicit def dcm[A, Raw] = new DCMultiplier[A, Raw] {}

  trait RDMultiplier[A, Raw] extends MatrixMultiplier[A, Raw, DenseRow, Dense, DenseRow] {
    def gemm(alpha: A, beta: A, m1: DenseRow[A, Raw], m2: Dense[A, Raw], ret: DenseRow[A, Raw]) =
      genericGemm(alpha, beta, m1, m2, ret)
  }
  implicit def rdm[A, Raw] = new RDMultiplier[A, Raw] {}
}


// Builders
trait DenseBuilders {
  
  implicit def dense[A, Raw](implicit sb: ScalarData.Builder[A, Raw]) = new MatrixBuilder[A, Raw, Dense] {
    def zeros(numRows: Int, numCols: Int) = {
      require(numRows > 0 && numCols > 0, "Cannot build matrix with non-positive dimensions (%d, %d)".format(numRows, numCols))
      val nr = numRows
      val nc = numCols
      new Dense[A, Raw] {
        val data: ScalarData[A, Raw] = sb.build(nr*nc)
        val scalar: ScalarOps[A] = data.scalar
        val numRows = nr
        val numCols = nc
      }
    }
  }
  
  implicit def denseRow[A, Raw](implicit sb: ScalarData.Builder[A, Raw]) = new MatrixBuilder[A, Raw, DenseRow] {
    def zeros(numRows: Int, numCols: Int) = {
      require(numRows > 0 && numCols > 0, "Cannot build matrix with non-positive dimensions (%d, %d)".format(numRows, numCols))
      require(numRows == 1, "Cannot build row matrix with %d rows".format(numRows))
      val nc = numCols
      new DenseRow[A, Raw] {
        val data: ScalarData[A, Raw] = sb.build(1*nc)
        val scalar: ScalarOps[A] = data.scalar
        val numRows = 1
        val numCols = nc
      }
    }
  }
  
  implicit def denseCol[A, Raw](implicit sb: ScalarData.Builder[A, Raw]) = new MatrixBuilder[A, Raw, DenseCol] {
    def zeros(numRows: Int, numCols: Int) = {
      require(numRows > 0 && numCols > 0, "Cannot build matrix with non-positive dimensions (%d, %d)".format(numRows, numCols))
      require(numCols == 1, "Cannot build column matrix with %d cols".format(numCols))
      val nr = numRows
      new DenseCol[A, Raw] {
        val data: ScalarData[A, Raw] = sb.build(nr*1)
        val scalar: ScalarOps[A] = data.scalar
        val numRows = nr
        val numCols = 1
      }
    }
  }
  
  class DenseBuilderExtras[A, Raw](implicit sb: ScalarData.Builder[A, Raw]) {
    def zeros(numRows: Int, numCols: Int): Dense[A, Raw] = {
      dense.zeros(numRows, numCols)
    }

    def tabulate(numRows: Int, numCols: Int)(f: (Int, Int) => A): Dense[A, Raw] = {
      val m = dense.zeros(numRows, numCols)
      m.tabulate(f)
      m
    }

    def eye(numRows: Int): Dense[A, Raw] = {
      val m = dense.zeros(numRows, numRows)
      m.tabulate { (i, j) => if (i == j) m.scalar.one else m.scalar.zero }
      m
    }
    
    def col(elems: A*): DenseCol[A, Raw] = {
      val m = denseCol.zeros(elems.size, 1)
      for (i <- 0 until m.numRows) m(i, 0) = elems(i)
      m
    }

    def row(elems: A*): DenseRow[A, Raw] = {
      val m = denseRow.zeros(1, elems.size)
      for (j <- 0 until m.numCols) m(0, j) = elems(j)
      m
    }

    def fromRows(row1: DenseRow[A, Raw], rows: DenseRow[A, Raw]*): Dense[A, Raw] = {
      require(row1.numRows == 1 && rows.forall(_.numRows == 1))
      require(rows.forall(_.numCols == row1.numCols))
      val ret = dense.zeros(1 + rows.size, row1.numCols)
      ret(0, ::) = row1
      for (i <- rows.indices)
        ret(i+1, ::) = rows(i)
      ret
    }
  }
  
  val denseRealFlt = new DenseBuilderExtras[Float, Float]
  val denseRealDbl = new DenseBuilderExtras[Double, Double]
  val denseComplexFlt = new DenseBuilderExtras[Complex, Float]
  val denseComplexDbl = new DenseBuilderExtras[Complex, Double]
}


