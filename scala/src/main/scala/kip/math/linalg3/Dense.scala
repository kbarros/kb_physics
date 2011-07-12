package kip.math.linalg3

import kip.math.Complex
import kip.math.linalg2.{ScalarOps, ScalarData}

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

trait DenseAdders1 {
  trait DDAdder[A, Raw] extends MatrixAdder[A, Raw, Dense, Dense, Dense] {
    def addInPlace(neg1: Boolean, m1: Dense[A, Raw], neg2: Boolean, m2: Dense[A, Raw], ret: Dense[A, Raw]) {
      println("dense dense adder")
    }
  }
  implicit def dda[A, Raw] = new DDAdder[A, Raw] {}
}

trait DenseAdders0 extends DenseAdders1 {
  trait DRAdder[A, Raw] extends MatrixAdder[A, Raw, Dense, DenseRow, DenseRow] {
    def addInPlace(neg1: Boolean, m1: Dense[A, Raw], neg2: Boolean, m2: DenseRow[A, Raw], ret: DenseRow[A, Raw]) {
    }
  }
  implicit def dra[A, Raw] = new DRAdder[A, Raw] {}

  trait RDAdder[A, Raw] extends MatrixAdder[A, Raw, DenseRow, Dense, DenseRow] {
    def addInPlace(neg1: Boolean, m1: DenseRow[A, Raw], neg2: Boolean, m2: Dense[A, Raw], ret: DenseRow[A, Raw]) {
    }
  }
  implicit def rda[A, Raw] = new RDAdder[A, Raw] {}
}

trait DenseAdders extends DenseAdders0{
  trait RRAdder[A, Raw] extends MatrixAdder[A, Raw, DenseRow, DenseRow, DenseRow] {
    override def addInPlace(neg1: Boolean, m1: DenseRow[A, Raw], neg2: Boolean, m2: DenseRow[A, Raw], ret: DenseRow[A, Raw]) {
      println("r r adder")
    }
  }
  implicit def rra[A, Raw] = new RRAdder[A, Raw] {}
}


// Builders
trait DenseBuilders {
  implicit def dense[A, Raw](implicit sb: ScalarData.Builder[A, Raw]) = new MatrixBuilder[A, Raw, Dense] {
    def zeros(numRows: Int, numCols: Int) = {
      val nr = numRows
      val nc = numCols
      new Dense[A, Raw] {
        val data: ScalarData[A, Raw] = sb.build(numRows*numCols)
        val scalar: ScalarOps[A] = data.scalar
        val numRows = nr
        val numCols = nc
      }
    }
  }
  
  implicit def denseRow[A, Raw](implicit sb: ScalarData.Builder[A, Raw]) = new MatrixBuilder[A, Raw, DenseRow] {
    def zeros(numRows: Int, numCols: Int) = {
      val nr = numRows
      val nc = numCols
      new DenseRow[A, Raw] {
        val data: ScalarData[A, Raw] = sb.build(numRows*numCols)
        val scalar: ScalarOps[A] = data.scalar
        val numRows = nr
        val numCols = nc
      }
    }
  }
  
}


