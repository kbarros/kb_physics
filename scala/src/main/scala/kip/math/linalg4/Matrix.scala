package kip.math.linalg4

import kip.math.Complex


// TODO: make GenMatrix that specializes S#A for relevant methods

trait Matrix[S <: Scalar, +Repr[S2 <: Scalar] <: Matrix[S2, Repr]] { self: Repr[S] =>
  val scalar: ScalarOps[S]
  
  def numRows: Int
  def numCols: Int
  def apply(i: Int, j: Int): S#A
  def update(i: Int, j: Int, x: S#A)
  def indices: Traversable[(Int, Int)]

  
  def checkKey(i: Int, j: Int) {
    require(0 <= i && i < numRows && 0 <= j && j < numCols, "Matrix indices out of bounds: [%d %d](%d %d)".format(numRows, numCols, i, j))
  }

  def transform(f: S#A => S#A): this.type = {
    indices.foreach { case(i, j) => this(i, j) = f(this(i, j)) }
    this
  }
  
  // The extra type parameter A2 and evidence parameter S2 are for type inference
  def map[A2, S2 <: Scalar{type A=A2}, That[S <: Scalar] >: Repr[S] <: Matrix[S, That]]
      (f: S#A => S2#A)(implicit ev: S2, mb: MatrixBuilder[S2, That]): That[S2] = {
    val ret = mb.zeros(numRows, numCols)
    indices.foreach { case(i, j) => ret(i, j) = f(this(i, j)) }
    ret
  }

  def duplicate[That[S <: Scalar] >: Repr[S] <: Matrix[S, That]](implicit mb: MatrixBuilder[S, That]): That[S] = {
    val ret = mb.zeros(numRows, numCols)
    indices.foreach { case(i, j) => ret(i, j) = this(i, j) }
    ret
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


  def *[Repr1[S <: Scalar] >: Repr[S], Repr2[S <: Scalar] <: Matrix[S, Repr2], Repr3[S <: Scalar] <: Matrix[S, Repr3]]
        (that: Repr2[S])
        (implicit mm: MatrixMultiplier[S, Repr1, Repr2, Repr3],
                  mb: MatrixBuilder[S, Repr3]): Repr3[S] = {
    require(numCols == that.numRows,
            "Can't multiply matrices of shape [%d, %d] * [%d, %d]".format(numRows, numCols, that.numRows, that.numCols))
    val ret = mb.zeros(numRows, that.numCols)
    mm.gemm(scalar.one, scalar.zero, this, that, ret)
    ret
  }

  def +[Repr1[S <: Scalar] >: Repr[S], Repr2[S <: Scalar] <: Matrix[S, Repr2], Repr3[S <: Scalar] <: Matrix[S, Repr3]]
        (that: Repr2[S])
        (implicit ma: MatrixAdder[S, Repr1, Repr2, Repr3],
                  mb: MatrixBuilder[S, Repr3]): Repr3[S] = {
    require(numRows == that.numRows && numCols == that.numCols,
            "Can't add matrices of shape [%d, %d] + [%d, %d]".format(numRows, numCols, that.numRows, that.numCols))
    val ret = mb.zeros(numRows, numCols)
    ma.addInPlace(false, this, that, ret)
    ret
  }

  def -[Repr1[S <: Scalar] >: Repr[S], Repr2[S <: Scalar] <: Matrix[S, Repr2], Repr3[S <: Scalar] <: Matrix[S, Repr3]]
        (that: Repr2[S])
        (implicit ma: MatrixAdder[S, Repr1, Repr2, Repr3],
                  mb: MatrixBuilder[S, Repr3]): Repr3[S] = {
    require(numRows == that.numRows && numCols == that.numCols,
            "Can't subtract matrices of shape [%d, %d] - [%d, %d]".format(numRows, numCols, that.numRows, that.numCols))
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

object MatrixAdder extends DenseAdders
trait MatrixAdder[S <: Scalar, Repr1[_ <: Scalar], Repr2[_ <: Scalar], Repr3[_ <: Scalar]] {
  def addInPlace(sub: Boolean, m1: Repr1[S], m2: Repr2[S], ret: Repr3[S])
}

object MatrixMultiplier extends DenseMultipliers
trait MatrixMultiplier[S <: Scalar, Repr1[_ <: Scalar], Repr2[_ <: Scalar], Repr3[_ <: Scalar]] {
  def gemm(alpha: S#A, beta: S#A, m1: Repr1[S], m2: Repr2[S], ret: Repr3[S])
}

object MatrixBuilder extends DenseBuilders
trait MatrixBuilder[S <: Scalar, Repr[_ <: Scalar]] {
  def zeros(numRows: Int, numCols: Int): Repr[S]
}


object Test extends App {

  val m1 = MatrixBuilder.denseRealDbl.zeros(4, 4)
  val m2 = MatrixBuilder.denseRealDbl.zeros(4, 4)
  val m3 = m1+m2
  println(m3)
  val m4 = MatrixBuilder.denseRealDbl.zeros(4, 4)
  m3.map(_.toDouble).map(_+Complex.I)

  
//  import kip.util.Util.time2
//  
//  time2("Eigenvalues") {
//    import MatrixBuilder.denseComplexDbl._
//    val n = 2000
//    val m3 = tabulate(n, n) { case (i, j) => i + 2*j }
//    val x = tabulate(n, 1) { case (i, j) => i }
//    val (v, w) = m3.eig
//    m3 * w(::,0) / v(0) - w(::,0)
//  }
}

