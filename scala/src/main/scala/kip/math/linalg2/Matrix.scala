package kip.math.linalg2

import kip.math.Complex



object MatrixImpl {

  object Builder {
    implicit def dense[@specialized(Float, Double) A, @specialized(Float, Double) Raw]
    (implicit sb: ScalarData.Builder[A, Raw]) = new Builder[A, Raw, Dense[A, Raw]] {
      def build(numRows: Int, numCols: Int) = new Dense(numRows, numCols)
    }
    
  }
  
  trait Builder[@specialized(Float, Double) A, @specialized(Float, Double) Raw, Impl <: MatrixImpl[A, Raw]] {
    def build(numRows: Int, numCols: Int): Impl
  }
  
  
  class Dense[@specialized(Float, Double) A, @specialized(Float, Double) Raw]
       (val numRows: Int, val numCols: Int)(implicit sb: ScalarData.Builder[A, Raw]) extends MatrixImpl[A, Raw] {
    val data: ScalarData[A, Raw] = sb.build(numRows*numCols)
    val scalar: ScalarOps[A] = data.scalar

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
    
    def transposeInPlace() {
      for (j <- 1 until numCols; i <- 0 until j) {
        val temp = this(i, j)
        this(i, j) = this(j, i)
        this(j, i) = temp
      }
    }
  }
}



trait MatrixImpl[@specialized(Float, Double) A, @specialized(Float, Double) Raw] {
  val scalar: ScalarOps[A]
  
  val numRows: Int
  val numCols: Int
  def apply(i: Int, j: Int): A
  def update(i: Int, j: Int, x: A)
  
  def indices: Iterator[(Int, Int)]
  
  def copyTo[That <: MatrixImpl[A, Raw]](that: That)(implicit ev: this.type <:< That) { // TODO: understand evidence parameter
    indices.foreach { case(i, j) => that(i, j) = this(i, j) }
  }

  def mapInPlace(f: A => A): Unit = indices.foreach { case(i, j) => this(i, j) = f(this(i, j)) }
  def multiplyInPlace(x: A): Unit = mapInPlace { scalar.mul(_, x) }
  def divideInPlace(x: A): Unit = mapInPlace { scalar.div(_, x) }
  def negateInPlace(): Unit = mapInPlace { scalar.neg(_) }
  def conjugateInPlace(): Unit = mapInPlace { scalar.conj(_) }
  def transposeInPlace(): Unit
}



// --------------
// Matrix Adder

trait MatrixAdder[@specialized(Float, Double) A, Raw, Impl1 <: MatrixImpl[A, Raw], Impl2 <: MatrixImpl[A, Raw], That] {
  def add(m1: Matrix[A, Raw, Impl1], m2: Matrix[A, Raw, Impl2]): That
}

// -----------

trait Matrix[@specialized(Float, Double) A, @specialized(Float, Double) Raw, Impl <: MatrixImpl[A, Raw]] {
  val impl: Impl
  
  def numRows: Int = impl.numRows
  def numCols: Int = impl.numCols
  def apply(i: Int, j: Int): A = impl.apply(i, j)
  def update(i: Int, j: Int, x: A): Unit = impl.update(i, j, x)
  
  def indices: Iterator[(Int, Int)] = impl.indices
  
  // def map[A2, Raw2, Impl2, That <: Matrix[A2, Raw2, Impl2]](f: A => A2)(implicit mb: MatrixBuilder[A2, Raw2, Impl, That]): That = {
  //   val ret = mb.zeros(numRows, numCols)
  //   indices.foreach { case (i, j) => ret(i, j) = f(this(i, j)) }
  //   ret
  // }
  
  def +[Impl2 <: MatrixImpl[A, Raw], That](that: Matrix[A, Raw, Impl2])(implicit ma: MatrixAdder[A, Raw, Impl, Impl2, That]) {
    
  }

  def clone[That <: Matrix[A, Raw, Impl]](implicit mb: MatrixBuilder[A, Raw, Impl, That]) = {
    val ret = mb.zeros(numRows, numCols)
    impl.copyTo(ret.impl)
    ret
  }


  def *[That <: Matrix[A, Raw, Impl]](x: A)(implicit mb: MatrixBuilder[A, Raw, Impl, That]) = {
    val ret = clone
    ret.impl.multiplyInPlace(x)
    ret
  }

  def /[That <: Matrix[A, Raw, Impl]](x: A)(implicit mb: MatrixBuilder[A, Raw, Impl, That]) = {
    val ret = clone
    ret.impl.divideInPlace(x)
    ret
  }

  def unary_-[That <: Matrix[A, Raw, Impl]](implicit mb: MatrixBuilder[A, Raw, Impl, That]) = {
    val ret = clone
    ret.impl.negateInPlace()
    ret
  }

  def tran[That <: Matrix[A, Raw, Impl]](implicit mb: MatrixBuilder[A, Raw, Impl, That]) = {
    val ret = clone
    ret.impl.transposeInPlace()
    ret
  }
  
  def conj[That <: Matrix[A, Raw, Impl]](implicit mb: MatrixBuilder[A, Raw, Impl, That]) = {
    val ret = clone
    ret.impl.conjugateInPlace()
    ret
  }
  
  def dag[That <: Matrix[A, Raw, Impl]](implicit mb: MatrixBuilder[A, Raw, Impl, That]) = {
    val ret = tran
    ret.impl.conjugateInPlace()
    ret
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


//  override def clone()(implicit b: MBuilder[A, Raw, Shape, Repr]: Repr = {
//    new DenseMatrix(numRows, numCols, scalar.clone(data))
//  }
/*  
  // TODO: Introduce "Repr" parameter to get most specific return type
  def map[A0, Raw0, That](f: A => A0)(implicit b: MBuilder[A0, Raw0, Shape, That]): That = {
    val ret = b.build
    

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
    val ret = matrix.zeros(numRows, numCols)
    DenseMatrix.QRSolve(ret, this, that, false)
    ret
  }
*/
}



object Matrix {
  object Dense {
    trait RealFlt    extends Matrix[Float,   Float,  MatrixImpl.Dense[Float,   Float]]
    trait RealDbl    extends Matrix[Double,  Double, MatrixImpl.Dense[Double,  Double]]
    trait ComplexFlt extends Matrix[Complex, Float,  MatrixImpl.Dense[Complex, Float]]
    trait ComplexDbl extends Matrix[Complex, Double, MatrixImpl.Dense[Complex, Double]]
    
  }
}


object MatrixBuilder {
  // // Generic builder -- use only specific builders instead
  // implicit def generic[@specialized(Float, Double) A, @specialized(Float, Double) Raw, Impl <: MatrixImpl[A, Raw]]
  // (implicit mib: MatrixImpl.Builder[A, Raw, Impl]) = new MatrixBuilder[A, Raw, Impl, Matrix[A, Raw, Impl]] {
  //   def zeros(numRows: Int, numCols: Int) = new { val impl = mib.build(numRows, numCols) } with Matrix[A, Raw, Impl] 
  // }

  implicit val DenseRealFlt = new Dense[Float, Float, MatrixImpl.Dense[Float, Float], Matrix.Dense.RealFlt] {
    def zeros(numRows: Int, numCols: Int) = new {
      val impl = new MatrixImpl.Dense[Float, Float](numRows, numCols)
    } with Matrix.Dense.RealFlt
  }
  implicit val DenseRealDbl = new Dense[Double, Double, MatrixImpl.Dense[Double, Double], Matrix.Dense.RealDbl] {
    def zeros(numRows: Int, numCols: Int) = new {
      val impl = new MatrixImpl.Dense[Double, Double](numRows, numCols)
    } with Matrix.Dense.RealDbl
  }
  implicit val DenseComplexFlt = new Dense[Complex, Float, MatrixImpl.Dense[Complex, Float], Matrix.Dense.ComplexFlt] {
    def zeros(numRows: Int, numCols: Int) = new {
      val impl = new MatrixImpl.Dense[Complex, Float](numRows, numCols)
    } with Matrix.Dense.ComplexFlt
  }
  implicit val DenseComplexDbl = new Dense[Complex, Double, MatrixImpl.Dense[Complex, Double], Matrix.Dense.ComplexDbl] {
    def zeros(numRows: Int, numCols: Int) = new {
      val impl = new MatrixImpl.Dense[Complex, Double](numRows, numCols)
    } with Matrix.Dense.ComplexDbl
  }
  
  trait Dense[@specialized(Float, Double) A, @specialized(Float, Double) Raw, Impl <: MatrixImpl[A, Raw], Repr <: Matrix[A, Raw, Impl]]
      extends MatrixBuilder[A, Raw, Impl, Repr] {
    def tabulate(numRows: Int, numCols: Int)(f: (Int, Int) => A): Repr = {
      val m = zeros(numRows, numCols)
      m.indices.foreach { case(i, j) => m(i, j) = f(i, j) }
      m
    }
  }
}

trait MatrixBuilder[@specialized(Float, Double) A, @specialized(Float, Double) Raw, Impl <: MatrixImpl[A, Raw], Repr <: Matrix[A, Raw, Impl]] {
  def zeros(numRows: Int, numCols: Int): Repr
}


object Test extends App {
  val mb1 = implicitly[MatrixBuilder[Double, Double, MatrixImpl.Dense[Double, Double], _]]
  println(mb1.getClass)

  def m1 = MatrixBuilder.DenseRealDbl.zeros(3, 3)
  println(m1)
  
  val m2 = m1.clone
  m2(1, 2) = 3
  println(m1*2)
  println(m2.dag.tran)
  
  import Complex._
  def m3 = MatrixBuilder.DenseComplexFlt.tabulate(3, 3){ case (i, j) => j+I }
  println(m3)
}
