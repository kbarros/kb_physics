package kip.math.linalg4

import kip.math.Complex


// TODO: make GenMatrix that specializes S#A for relevant methods

trait Matrix[S <: Scalar, +Repr[S2 <: Scalar] <: Matrix[S2, Repr]] { self: Repr[S] =>
  val scalar: ScalarOps[S]
  
  def numRows: Int
  def numCols: Int
  def apply(i: Int, j: Int): S#A
  def update(i: Int, j: Int, x: S#A)
  def indices: Iterator[(Int, Int)]

  def checkKey(i: Int, j: Int) {
    require(0 <= i && i < numRows && 0 <= j && j < numCols, "Matrix indices out of bounds: [%d %d](%d %d)".format(numRows, numCols, i, j))
  }

  def mapTo[S2 <: Scalar, That[S <: Scalar] >: Repr[S] <: Matrix[S, That]](f: S#A => S2#A, that: That[S2]): Unit =
    indices.foreach { case(i, j) => that(i, j) = f(this(i, j)) }
  def copyTo[That[S <: Scalar] >: Repr[S] <: Matrix[S, That]](that: That[S]): Unit =
    mapTo(identity, that)
  def mulTo[That[S <: Scalar] >: Repr[S] <: Matrix[S, That]](x: S#A, that: That[S]): Unit =
    mapTo(scalar.mul(_, x), that)
  def divTo[That[S <: Scalar] >: Repr[S] <: Matrix[S, That]](x: S#A, that: That[S]): Unit =
    mapTo(scalar.div(_, x), that)
  def negTo[That[S <: Scalar] >: Repr[S] <: Matrix[S, That]](that: That[S]): Unit =
    mapTo(scalar.neg(_), that)
  def conjTo[That[S <: Scalar] >: Repr[S] <: Matrix[S, That]](that: That[S]): Unit =
    mapTo(scalar.conj(_), that)

  def map[S2 <: Scalar, That[S <: Scalar] >: Repr[S] <: Matrix[S, That]]
      (f: S#A => S2#A)(implicit mb: MatrixBuilder[S2, That]): That[S2] = {
    val ret = mb.zeros(numRows, numCols)
    mapTo(f, ret)
    ret
  }

  def duplicate[That[S <: Scalar] >: Repr[S] <: Matrix[S, That]](implicit mb: MatrixBuilder[S, That]): That[S] = {
    val ret = mb.zeros(numRows, numCols)
    copyTo(ret)
    ret
  }
  
  def *[That[S <: Scalar] >: Repr[S] <: Matrix[S, That]](x: S#A)(implicit mb: MatrixBuilder[S, That]): That[S] = {
    val ret = duplicate(mb)
    ret.mulTo(x, ret)
    ret
  }

  def /[That[S <: Scalar] >: Repr[S] <: Matrix[S, That]](x: S#A)(implicit mb: MatrixBuilder[S, That]): That[S] = {
    val ret = duplicate(mb)
    ret.divTo(x, ret)
    ret
  }

  def unary_-[That[S <: Scalar] >: Repr[S] <: Matrix[S, That]](implicit mb: MatrixBuilder[S, That]): That[S] = {
    val ret = duplicate(mb)
    ret.negTo(ret)
    ret
  }

  def conj[That[S <: Scalar] >: Repr[S] <: Matrix[S, That]](implicit mb: MatrixBuilder[S, That]): That[S] = {
    val ret = duplicate(mb)
    ret.conjTo(ret)
    ret
  }


  def *[Repr1[S <: Scalar] >: Repr[S], Repr2[S <: Scalar] <: Matrix[S, Repr2], Repr3[S <: Scalar] <: Matrix[S, Repr3]]
        (that: Repr2[S])
        (implicit mm: MatrixMultiplier[S, Repr1, Repr2, Repr3],
                  mb: MatrixBuilder[S, Repr3]): Repr3[S] = {
    require(numCols == that.numRows,
            "Can't multiply matrices of dimensions (%d, %d) and (%d, %d)".format(numRows, numCols, that.numRows, that.numCols))
    val ret = mb.zeros(numRows, that.numCols)
    mm.gemm(scalar.one, scalar.zero, this, that, ret)
    ret
  }

  def +[Repr1[S <: Scalar] >: Repr[S], Repr2[S <: Scalar] <: Matrix[S, Repr2], Repr3[S <: Scalar] <: Matrix[S, Repr3]]
        (that: Repr2[S])
        (implicit ma: MatrixAdder[S, Repr1, Repr2, Repr3],
                  mb: MatrixBuilder[S, Repr3]): Repr3[S] = {
    require(numRows == that.numRows && numCols == that.numCols,
            "Can't add matrices of dimensions (%d, %d) and (%d, %d)".format(numRows, numCols, that.numRows, that.numCols))
    val ret = mb.zeros(numRows, numCols)
    ma.addInPlace(false, this, that, ret)
    ret
  }

  def -[Repr1[S <: Scalar] >: Repr[S], Repr2[S <: Scalar] <: Matrix[S, Repr2], Repr3[S <: Scalar] <: Matrix[S, Repr3]]
        (that: Repr2[S])
        (implicit ma: MatrixAdder[S, Repr1, Repr2, Repr3],
                  mb: MatrixBuilder[S, Repr3]): Repr3[S] = {
    require(numRows == that.numRows && numCols == that.numCols,
            "Can't subtract matrices of dimensions (%d, %d) and (%d, %d)".format(numRows, numCols, that.numRows, that.numCols))
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
  
  // TODO: Make this more reasonable
  m3.mapTo[Scalar.RealDbl, Dense](_.toDouble, m4)

//  println(implicitly[MatrixAdder[Float, Float, DenseRow, DenseRow, Dense]])
}



/*
// why compiles?
trait Matrix[-Repr <: Matrix[Repr]]
trait Dense extends Matrix[Dense]
trait DenseRow extends Matrix[DenseRow] with Dense
*/
