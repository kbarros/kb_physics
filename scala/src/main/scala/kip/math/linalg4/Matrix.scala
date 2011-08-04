package kip.math.linalg4


object MatrixDims {
  def checkDims(numRows: Int, numCols: Int) {
    require(numRows > 0 && numCols > 0,
        "Cannot build matrix with non-positive dimensions [%d, %d]".format(numRows, numCols))
  }
  
  def checkKey(m: MatrixDims, i: Int, j: Int) {
    require(0 <= i && i < m.numRows && 0 <= j && j < m.numCols,
        "Matrix indices out of bounds: [%d %d](%d %d)".format(m.numRows, m.numCols, i, j))
  }
  
  def checkAdd(m1: MatrixDims, m2: MatrixDims) {
    require(m1.numRows == m2.numRows && m1.numCols == m2.numCols,
        "Cannot add/subtract matrices of shape [%d, %d] +- [%d, %d]".format(m1.numRows, m1.numCols, m2.numRows, m2.numCols))
  }
  def checkAddTo(m1: MatrixDims, m2: MatrixDims, ret: MatrixDims) {
    checkAdd(m1, m2)
    require(
        ret.numRows == m1.numRows &&
        ret.numRows == m2.numRows &&
        ret.numCols == m1.numCols &&
        ret.numCols == m2.numCols,
        "Cannot add/subtract matrices of shape: [%d, %d] +- [%d, %d] -> [%d, %d]".format(
            m1.numRows, m1.numCols, m2.numRows, m2.numCols, ret.numRows, ret.numCols))
  }

  def checkMul(m1: MatrixDims, m2: MatrixDims) {
    require(m1.numCols == m2.numRows,
            "Cannot multiply matrices of shape [%d, %d] * [%d, %d]".format(m1.numRows, m1.numCols, m2.numRows, m2.numCols))
  }
  def checkMulTo(m1: MatrixDims, m2: MatrixDims, ret: MatrixDims) {
    checkMul(m1, m2)
    require(
        ret.numRows == m1.numRows &&
        m1.numCols == m2.numRows &&
        m2.numCols == ret.numCols, "Cannot multiply matrices of shape: [%d, %d] * [%d, %d] -> [%d, %d]".format(
            m1.numRows, m1.numCols, m2.numRows, m2.numCols, ret.numRows, ret.numCols))
  }
  
  def checkDot(m1: MatrixDims, m2: MatrixDims) {
    checkMul(m1, m2)
    require(m1.numRows == 1 && m2.numCols == 1,
            "Dot product expects row and column vectors, found [%d, %d] * [%d, %d]".format(m1.numRows, m1.numCols, m2.numRows, m2.numCols))
  }
}


trait MatrixDims {
  def numRows: Int
  def numCols: Int
}


// TODO: Specialize S#A for apply/update methods

trait Matrix[S <: Scalar, +Repr[S2 <: Scalar] <: Matrix[S2, Repr]] extends MatrixDims { self: Repr[S] =>
  val scalar: ScalarOps[S]
  
  def apply(i: Int, j: Int): S#A
  def update(i: Int, j: Int, x: S#A)
  def transform(f: S#A => S#A): this.type
  def dispose() {}
  
  def duplicate[That[S <: Scalar] >: Repr[S] <: Matrix[S, That]](implicit mb: MatrixBuilder[S, That]): That[S] = {
    mb.duplicate(this)
  }
  
  // The parameter A2 is for type inference only
  def map[A2, S2 <: Scalar{type A=A2}, That[T <: Scalar] >: Repr[T] <: Matrix[T, That]]
      (f: S#A => S2#A)(implicit scalar2: ScalarOps[S2], mb: MatrixBuilder[S2, That]): That[S2] = {
    mb.map(this)(f)
  }
  
  def *[That[S <: Scalar] >: Repr[S] <: Matrix[S, That]](x: S#A)(implicit mb: MatrixBuilder[S, That]): That[S] = {
    duplicate(mb).transform(scalar.mul(_, x))
  }

  def /[That[S <: Scalar] >: Repr[S] <: Matrix[S, That]](x: S#A)(implicit mb: MatrixBuilder[S, That]): That[S] = {
    duplicate(mb).transform(scalar.div(_, x))
  }

  def unary_-[That[S <: Scalar] >: Repr[S] <: Matrix[S, That]](implicit mb: MatrixBuilder[S, That]): That[S] = {
    duplicate(mb).transform(scalar.neg(_))
  }

  def conj[That[S <: Scalar] >: Repr[S] <: Matrix[S, That]](implicit mb: MatrixBuilder[S, That]): That[S] = {
    duplicate(mb).transform(scalar.conj(_))
  }

  def tran[That[S <: Scalar] >: Repr[S] <: Matrix[S, That]](implicit mb: MatrixBuilder[S, That]): That[S] = {
    mb.transpose(this)
  }
  
  def dag[That[S <: Scalar] >: Repr[S] <: Matrix[S, That]](implicit mb: MatrixBuilder[S, That]): That[S] = {
    tran(mb).transform(scalar.conj(_))
  }

  def dot[Repr1[S <: Scalar] >: Repr[S], Repr2[S <: Scalar] <: Matrix[S, Repr2], Repr3[S <: Scalar] <: Matrix[S, Repr3]]
        (that: Repr2[S])
        (implicit mm: MatrixMultiplier[S, Repr1, Repr2, Repr3],
                  mb: MatrixBuilder[S, Repr3]): S#A = {
    MatrixDims.checkDot(this, that)
    val m3 = mb.zeros(1, 1)
    mm.gemm(scalar.one, scalar.zero, this, that, m3)
    val ret = m3(0, 0)
    m3.dispose()
    ret
  }

  def *[Repr1[S <: Scalar] >: Repr[S], Repr2[S <: Scalar] <: Matrix[S, Repr2], Repr3[S <: Scalar] <: Matrix[S, Repr3]]
        (that: Repr2[S])
        (implicit mm: MatrixMultiplier[S, Repr1, Repr2, Repr3],
                  mb: MatrixBuilder[S, Repr3]): Repr3[S] = {
    MatrixDims.checkMul(this, that)
    val ret = mb.zeros(numRows, that.numCols)
    mm.gemm(scalar.one, scalar.zero, this, that, ret)
    ret
  }

  def +[Repr1[S <: Scalar] >: Repr[S], Repr2[S <: Scalar] <: Matrix[S, Repr2], Repr3[S <: Scalar] <: Matrix[S, Repr3]]
        (that: Repr2[S])
        (implicit ma: MatrixAdder[S, Repr1, Repr2, Repr3],
                  mb: MatrixBuilder[S, Repr3]): Repr3[S] = {
    MatrixDims.checkAdd(this, that)
    val ret = mb.zeros(numRows, numCols)
    ma.addInPlace(false, this, that, ret)
    ret
  }

  def -[Repr1[S <: Scalar] >: Repr[S], Repr2[S <: Scalar] <: Matrix[S, Repr2], Repr3[S <: Scalar] <: Matrix[S, Repr3]]
        (that: Repr2[S])
        (implicit ma: MatrixAdder[S, Repr1, Repr2, Repr3],
                  mb: MatrixBuilder[S, Repr3]): Repr3[S] = {
    MatrixDims.checkAdd(this, that)
    val ret = mb.zeros(numRows, numCols)
    ma.addInPlace(true, this, that, ret)
    ret
  }
  
  override def toString = {
    val sb = new StringBuilder()
    val maxRows = 6
    val maxCols = 6
    val elemSpacing = 2
    
    var elemWidth = 8
    for (i <- 0 until math.min(numRows, maxRows)) {
      for (j <- 0 until math.min(numCols, maxCols)) {
        elemWidth = math.max(elemWidth, this(i, j).toString.size+elemSpacing)
      }
    }
    
    def writeStr(str: String) {
      val spaces = elemWidth-str.size
      val margin1 = spaces / 2
      val margin2 = spaces - margin1
      sb.append(Seq.fill(margin1)(' ').mkString)
      sb.append(str)
      sb.append(Seq.fill(margin2)(' ').mkString)
    }
    for (i <- 0 until math.min(numRows, maxRows)) {
      for (j <- 0 until math.min(numCols, maxCols)) {
        writeStr(this(i, j).toString)
      }
      if (i == 0 && numCols > maxCols)
        sb.append(" ... (%d Cols)".format(numCols))
      sb.append("\n")
    }
    if (numRows > maxRows) {
      writeStr(":")
      sb.append("\n")
      writeStr("(%d Rows)".format(numRows))
    }
    sb.toString
  }
}


trait MatrixBuilder[S <: Scalar, Repr[_ <: Scalar]] {
  def zeros(numRows: Int, numCols: Int): Repr[S]
  def duplicate(m: Repr[S]): Repr[S]
  def transpose(m: Repr[S]): Repr[S]
  def map[S0 <: Scalar](m: Repr[S0])(f: S0#A => S#A): Repr[S]
}

trait MatrixAdder[S <: Scalar, Repr1[_ <: Scalar], Repr2[_ <: Scalar], Repr3[_ <: Scalar]] {
  def addInPlace(sub: Boolean, m1: Repr1[S], m2: Repr2[S], ret: Repr3[S])
}

trait MatrixMultiplier[S <: Scalar, Repr1[_ <: Scalar], Repr2[_ <: Scalar], Repr3[_ <: Scalar]] {
  def gemm(alpha: S#A, beta: S#A, m1: Repr1[S], m2: Repr2[S], ret: Repr3[S])
}


object Test extends App {
  import Constructors.complexDbl._
  val m = dense(4, 4)
  m(0, 0) = 1

  println((m+m).map(_.re))
  
  import kip.util.Util.time2
  
  time2("Eigenvalues") {
    val n = 500
    val m3 = tabulate(n, n) { case (i, j) => i + 2*j }
    val m4 = m3 + m3.tran
    val x = tabulate(n, 1) { case (i, j) => i }
    val (v, w) = m3.eig
    println(m3 * w(::,0) / v(0) - w(::,0))
  }
}

