package kip.math.linalg3

import kip.math.Complex
import kip.math.linalg2.{ScalarOps, ScalarData}


// TODO: put implicits in object Dense 



// Matrix types

trait Dense[A, Raw] extends Matrix[A, Raw, Dense] {
  val data: ScalarData[A, Raw]
  
  def checkKey(i: Int, j: Int) {
    require(0 <= i && i < numRows)
    require(0 <= j && j < numCols)
  }

  def index(i: Int, j: Int) = {
    checkKey(i, j)
    i + j*numRows // fortran column major convention
  }
  
  def apply(i: Int, j: Int): A = data(index(i, j))
  def update(i: Int, j: Int, x: A) { data(index(i, j)) = x }
  
  def indices = (for (j <- 0 until numCols; i <- 0 until numRows) yield (i, j)).toIterator
}

trait DenseRow[A, Raw] extends Matrix[A, Raw, DenseRow] with Dense[A, Raw]
trait DenseCol[A, Raw] extends Matrix[A, Raw, DenseCol] with Dense[A, Raw]


// Adders

trait DenseAdders {
  private def genericAddInPlace[A, Raw](sub: Boolean, m1: Dense[A, Raw], m2: Dense[A, Raw], ret: Dense[A, Raw]) {
    require(ret.numRows == m1.numRows && ret.numRows == m2.numRows, "Mismatched rows: %d, %d, %d".format(m1.numRows, m2.numRows, ret.numRows))
    require(ret.numCols == m1.numCols && ret.numCols == m2.numCols, "Mismatched cols: %d, %d, %d".format(m1.numCols, m2.numCols, ret.numCols))
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

  trait DRMultiplier[A, Raw] extends MatrixMultiplier[A, Raw, Dense, DenseRow, DenseRow] {
    def gemm(alpha: A, beta: A, m1: Dense[A, Raw], m2: DenseRow[A, Raw], ret: DenseRow[A, Raw]) =
      genericGemm(alpha, beta, m1, m2, ret)
  }
  implicit def drm[A, Raw] = new DDMultiplier[A, Raw] {}
}


// Builders
trait DenseBuilders {

  implicit def dense[A, Raw](implicit sb: ScalarData.Builder[A, Raw]) = new MatrixBuilder[A, Raw, Dense] {
    def zeros(numRows: Int, numCols: Int) = {
      require(numRows > 0 && numCols > 0, "Can't build matrix with non-positive dimensions (%d, %d)".format(numRows, numCols))
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
      require(numRows > 0 && numCols > 0, "Can't build matrix with non-positive dimensions (%d, %d)".format(numRows, numCols))
      require(numRows == 1, "Can't build row matrix with %d rows".format(numRows))
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
      require(numRows > 0 && numCols > 0, "Can't build matrix with non-positive dimensions (%d, %d)".format(numRows, numCols))
      require(numCols == 1, "Can't build column matrix with %d rows".format(numCols))
      val nr = numRows
      new DenseCol[A, Raw] {
        val data: ScalarData[A, Raw] = sb.build(nr*1)
        val scalar: ScalarOps[A] = data.scalar
        val numRows = nr
        val numCols = 1
      }
    }
  }
  
}


