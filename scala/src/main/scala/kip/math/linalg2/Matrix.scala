package kip.math.linalg2

import kip.math.Complex


// ----------

object ScalarOps {
  implicit object RealDbl extends RealDbl
  implicit object RealFlt extends RealFlt
  implicit object ComplexDbl extends ComplexDbl

  trait RealDbl extends ScalarOps[Double] {
    def add(a: Double, b: Double): Double = a + b
    def sub(a: Double, b: Double): Double = a - b
    def mul(a: Double, b: Double): Double = a * b
    def div(a: Double, b: Double): Double = a / b
    def neg(a: Double): Double = -a
    def conj(a: Double): Double = a
    def zero: Double = 0.0
    def one: Double = 1.0
  }

  trait RealFlt extends ScalarOps[Float] {
    def add(a: Float, b: Float): Float = a + b
    def sub(a: Float, b: Float): Float = a - b
    def mul(a: Float, b: Float): Float = a * b
    def div(a: Float, b: Float): Float = a / b
    def neg(a: Float): Float = -a
    def conj(a: Float): Float = a
    def zero: Float = 0.0f
    def one: Float = 1.0f
  }

  trait ComplexDbl extends ScalarOps[Complex] {
    def add(a: Complex, b: Complex): Complex = a + b
    def sub(a: Complex, b: Complex): Complex = a - b
    def mul(a: Complex, b: Complex): Complex = a * b
    def div(a: Complex, b: Complex): Complex = a / b
    def neg(a: Complex): Complex = -a
    def conj(a: Complex): Complex = a.conj
    def zero: Complex = 0
    def one: Complex = 1
  }

}

trait ScalarOps[@specialized(Float, Double) T] {
  def add(a: T, b: T): T
  def sub(a: T, b: T): T
  def mul(a: T, b: T): T
  def div(a: T, b: T): T
  def neg(a: T): T
  def conj(a: T): T
  def zero: T
  def one: T
}


// -----------


object ScalarData {
  // Builder
  
  object Builder {
    implicit val RealFlt = new Builder[Float, Float] {
      def build(size: Int) = new RealFlt(size)
    }
    implicit val RealDbl = new Builder[Double, Double] {
      def build(size: Int) = new RealDbl(size)
    }
    implicit val ComplexFlt = new Builder[Complex, Float] {
      def build(size: Int) = new ComplexFlt(size)
    }
    implicit val ComplexDbl = new Builder[Complex, Double] {
      def build(size: Int) = new ComplexDbl(size)
    }
  }
  
  trait Builder[@specialized(Float, Double) A, @specialized(Float, Double) Raw] {
    def build(size: Int): ScalarData[A, Raw]
  }  
  
  // ScalarData Implementations                                          

  class RealFlt(size: Int) extends ScalarData[Float, Float] {
    val scalar = ScalarOps.RealFlt
    
    val raw = new Array[Float](size)
    def rawApply(i: Int): Float = raw(i)
    def rawUpdate(i: Int, x: Float) { raw(i) = x }
    
    def apply(i: Int): Float = raw(i)
    def update(i: Int, x: Float) = raw(i) = x
    override def madd(i0: Int, a1: ScalarData[Float, Float], i1: Int, a2: ScalarData[Float, Float], i2: Int) {
      raw(i0) += a1.rawApply(i1)*a2.rawApply(i2)
    }
  }

  class RealDbl(size: Int) extends ScalarData[Double, Double] {
    val scalar = ScalarOps.RealDbl
    
    val raw = new Array[Double](size)
    def rawApply(i: Int): Double = raw(i)
    def rawUpdate(i: Int, x: Double) { raw(i) = x }
    
    def apply(i: Int): Double = raw(i)
    def update(i: Int, x: Double) = raw(i) = x
    override def madd(i0: Int, a1: ScalarData[Double, Double], i1: Int, a2: ScalarData[Double, Double], i2: Int) {
      raw(i0) += a1.rawApply(i1)*a2.rawApply(i2)
    }
  }

  class ComplexFlt(size: Int) extends ScalarData[Complex, Float] {
    val scalar = ScalarOps.ComplexDbl
    
    val raw = new Array[Float](2*size)
    def rawApply(i: Int): Float = raw(i)
    def rawUpdate(i: Int, x: Float) { raw(i) = x }
    
    def apply(i: Int) = Complex(raw(2*i+0), raw(2*i+1))
    def update(i: Int, x: Complex) {
      raw(2*i+0) = x.re.toFloat
      raw(2*i+1) = x.im.toFloat
    }
    override def madd(i0: Int, a1: ScalarData[Complex, Float], i1: Int, a2: ScalarData[Complex, Float], i2: Int) {
      val x1_re = a1.rawApply(2*i1+0)
      val x1_im = a1.rawApply(2*i1+1)
      val x2_re = a2.rawApply(2*i2+0)
      val x2_im = a2.rawApply(2*i2+1)
      raw(2*i0+0) += x1_re*x2_re - x1_im*x2_im
      raw(2*i0+1) += x1_re*x2_im + x1_im*x2_re
    }
  }

  class ComplexDbl(size: Int) extends ScalarData[Complex, Double] {
    val scalar = ScalarOps.ComplexDbl
    
    val raw = new Array[Double](2*size)
    def rawApply(i: Int): Double = raw(i)
    def rawUpdate(i: Int, x: Double) { raw(i) = x }
    
    def apply(i: Int) = Complex(raw(2*i+0), raw(2*i+1))
    def update(i: Int, x: Complex) {
      raw(2*i+0) = x.re.toDouble
      raw(2*i+1) = x.im.toDouble
    }
    override def madd(i0: Int, a1: ScalarData[Complex, Double], i1: Int, a2: ScalarData[Complex, Double], i2: Int) {
      val x1_re = a1.rawApply(2*i1+0)
      val x1_im = a1.rawApply(2*i1+1)
      val x2_re = a2.rawApply(2*i2+0)
      val x2_im = a2.rawApply(2*i2+1)
      raw(2*i0+0) += x1_re*x2_re - x1_im*x2_im
      raw(2*i0+1) += x1_re*x2_im + x1_im*x2_re
    }
  }

}


abstract class ScalarData[@specialized(Float, Double) A, @specialized(Float, Double) Raw] {
  val scalar: ScalarOps[A]
  def rawApply(i: Int): Raw
  def rawUpdate(i: Int, x: Raw)
  def apply(i: Int): A
  def update(i: Int, x: A)
  def madd(i0: Int, a1: ScalarData[A, Raw], i1: Int, a2: ScalarData[A, Raw], i2: Int) {
    this(i0) = scalar.add(this(i0), scalar.mul(a1(i1), a2(i2)))
  }
}


// -----------


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
    
//    def map[A0, Raw0, That](f: A => A0)(implicit mb: MatrixImpl.Builder[A0, Raw0, That]): That
  }
}



trait MatrixImpl[@specialized(Float, Double) A, @specialized(Float, Double) Raw] {
  val numRows: Int
  val numCols: Int
  def apply(i: Int, j: Int): A
  def update(i: Int, j: Int, x: A)
  
  def indices: Iterator[(Int, Int)]
  def tabulate(f: (Int, Int) => A): Unit = indices.map { case(i, j) => this(i, j) = f(i, j) }
  def copyTo[That <: MatrixImpl[A, Raw]](that: That)(implicit ev: this.type <:< That) { // TODO: understand evidence parameter
    indices.map { case(i, j) => that(i, j) = this(i, j) }
  }
//  def map[A0, Raw0, That](f: A => A0)(implicit mb: MatrixImpl.Builder[A0, Raw0, That]): That
}



// --------------
// Matrix Adder

trait MatrixAdder[@specialized(Float, Double) A, Raw, Impl1 <: MatrixImpl[A, Raw], Impl2 <: MatrixImpl[A, Raw], That] {
  def add(m1: Matrix[A, Raw, Impl1], m2: Matrix[A, Raw, Impl2]): That
}

// -----------

object Matrix {
  trait RealFltDense    extends Matrix[Float,   Float,  MatrixImpl.Dense[Float,   Float]]
  trait RealDblDense    extends Matrix[Double,  Double, MatrixImpl.Dense[Double,  Double]]
  trait ComplexFltDense extends Matrix[Complex, Float,  MatrixImpl.Dense[Complex, Float]]
  trait ComplexDblDense extends Matrix[Complex, Double, MatrixImpl.Dense[Complex, Double]]
}

trait Matrix[@specialized(Float, Double) A, @specialized(Float, Double) Raw, Impl <: MatrixImpl[A, Raw]] {
  val impl: Impl
  
  def numRows: Int = impl.numRows
  def numCols: Int = impl.numCols
  def apply(i: Int, j: Int): A = impl.apply(i, j)
  def update(i: Int, j: Int, x: A): Unit = impl.update(i, j, x)
  
//  def map[A0, Raw0, That](f: A => A0)(implicit mb: MatrixImpl.Builder[A0, Raw0, That]): That = impl.map[A0, Raw0, That](f)(mb)
  
  def +[Impl2 <: MatrixImpl[A, Raw], That](that: Matrix[A, Raw, Impl2])(implicit ma: MatrixAdder[A, Raw, Impl, Impl2, That]) {
    
  }

  def clone[That <: Matrix[A, Raw, Impl]](implicit mb: MatrixBuilder[A, Raw, Impl, That]) = {
    val ret = mb.zeros(numRows, numCols)
    impl.copyTo(ret.impl)
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
    

.tabulate(numRows, numCols) { (i, j) => f(this(i, j)) }
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

*/
}



object MatrixBuilder {
  
  
  implicit val RealDblDense = new MatrixBuilder[Double, Double, MatrixImpl.Dense[Double, Double], Matrix.RealDblDense] {
    def zeros(numRows: Int, numCols: Int) = new { val impl = new MatrixImpl.Dense[Double, Double](numRows, numCols) } with Matrix.RealDblDense
  }

//  implicit def generic[@specialized(Float, Double) A, @specialized(Float, Double) Raw, Impl <: MatrixImpl[A, Raw]]
//  (implicit mib: MatrixImpl.Builder[A, Raw, Impl]) = new MatrixBuilder[A, Raw, Impl, Matrix[A, Raw, Impl]] {
//    def zeros(numRows: Int, numCols: Int) = new { val impl = mib.build(numRows, numCols) } with Matrix[A, Raw, Impl] 
//  }
}

trait MatrixBuilder[@specialized(Float, Double) A, @specialized(Float, Double) Raw, Impl <: MatrixImpl[A, Raw], Repr <: Matrix[A, Raw, Impl]] {
  def zeros(numRows: Int, numCols: Int): Repr
  
  def tabulate(numRows: Int, numCols: Int)(f: (Int, Int) => A): Repr = {
    val m = zeros(numRows, numCols)
    m.impl.tabulate(f)
    m
  }
}


object Test {
  def m = MatrixBuilder.RealDblDense.zeros(3, 3)
//  val mb = implicitly[MatrixBuilder[Double, Double, MatrixImpl.Dense[Double, Double], Matrix.RealDblDense]]
  val m2 = m.clone
}
