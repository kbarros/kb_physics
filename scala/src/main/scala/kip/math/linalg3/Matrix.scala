package kip.math.linalg3

import kip.math.Complex


trait Matrix[A, Raw, +Repr[C, D] <: Matrix[C, D, Repr]] { self: Repr[A, Raw] =>
  val scalar: ScalarOps[A]
  
  def numRows: Int
  def numCols: Int
  def apply(i: Int, j: Int): A
  def update(i: Int, j: Int, x: A)
  def indices: Iterator[(Int, Int)]

  def checkKey(i: Int, j: Int) {
    require(0 <= i && i < numRows && 0 <= j && j < numCols, "Matrix indices out of bounds: [%d %d](%d %d)".format(numRows, numCols, i, j))
  }

  def mapTo[A2, Raw2, That[C,D] >: Repr[C,D] <: Matrix[C, D, That]](f: A => A2, that: That[A2, Raw2]): Unit =
    indices.foreach { case(i, j) => that(i, j) = f(this(i, j)) }
  def copyTo[That[C,D] >: Repr[C,D] <: Matrix[C, D, That]](that: That[A, Raw]): Unit =
    mapTo(identity, that)
  def mulTo[That[C,D] >: Repr[C,D] <: Matrix[C, D, That]](x: A, that: That[A, Raw]): Unit =
    mapTo(scalar.mul(_, x), that)
  def divTo[That[C,D] >: Repr[C,D] <: Matrix[C, D, That]](x: A, that: That[A, Raw]): Unit =
    mapTo(scalar.div(_, x), that)
  def negTo[That[C,D] >: Repr[C,D] <: Matrix[C, D, That]](that: That[A, Raw]): Unit =
    mapTo(scalar.neg(_), that)
  def conjTo[That[C,D] >: Repr[C,D] <: Matrix[C, D, That]](that: That[A, Raw]): Unit =
    mapTo(scalar.conj(_), that)


  def map[A2, Raw2, That[C,D] >: Repr[C,D] <: Matrix[C, D, That]]
      (f: A => A2)(implicit mb: MatrixBuilder[A2, Raw2, That]): That[A2, Raw2] = {
    val ret = mb.zeros(numRows, numCols)
    mapTo(f, ret)
    ret
  }

  def duplicate[That[C,D] >: Repr[C,D] <: Matrix[C, D, That]](implicit mb: MatrixBuilder[A, Raw, That]): That[A, Raw] = {
    val ret = mb.zeros(numRows, numCols)
    copyTo(ret)
    ret
  }
  
  def *[That[C, D] >: Repr[C, D] <: Matrix[C, D, That]](x: A)(implicit mb: MatrixBuilder[A, Raw, That]): That[A, Raw] = {
    val ret = duplicate(mb)
    ret.mulTo(x, ret)
    ret
  }

  def /[That[C, D] >: Repr[C, D] <: Matrix[C, D, That]](x: A)(implicit mb: MatrixBuilder[A, Raw, That]): That[A, Raw] = {
    val ret = duplicate(mb)
    ret.divTo(x, ret)
    ret
  }

  def unary_-[That[C, D] >: Repr[C, D] <: Matrix[C, D, That]](implicit mb: MatrixBuilder[A, Raw, That]): That[A, Raw] = {
    val ret = duplicate(mb)
    ret.negTo(ret)
    ret
  }

  def conj[That[C, D] >: Repr[C, D] <: Matrix[C, D, That]](implicit mb: MatrixBuilder[A, Raw, That]): That[A, Raw] = {
    val ret = duplicate(mb)
    ret.conjTo(ret)
    ret
  }


  def *[Repr1[C, D] >: Repr[C, D], Repr2[C, D] <: Matrix[C, D, Repr2], Repr3[C, D] <: Matrix[C, D, Repr3]]
        (that: Repr2[A, Raw])
        (implicit mm: MatrixMultiplier[A, Raw, Repr1, Repr2, Repr3],
                  mb: MatrixBuilder[A, Raw, Repr3]): Repr3[A, Raw] = {
    require(numCols == that.numRows,
            "Can't multiply matrices of dimensions (%d, %d) and (%d, %d)".format(numRows, numCols, that.numRows, that.numCols))
    val ret = mb.zeros(numRows, that.numCols)
    mm.gemm(scalar.one, scalar.zero, this, that, ret)
    ret
  }

  def +[Repr1[C, D] >: Repr[C, D], Repr2[C, D] <: Matrix[C, D, Repr2], Repr3[C, D] <: Matrix[C, D, Repr3]]
        (that: Repr2[A, Raw])
        (implicit mm: MatrixAdder[A, Raw, Repr1, Repr2, Repr3],
                  mb: MatrixBuilder[A, Raw, Repr3]): Repr3[A, Raw] = {
    require(numRows == that.numRows && numCols == that.numCols,
            "Can't add matrices of dimensions (%d, %d) and (%d, %d)".format(numRows, numCols, that.numRows, that.numCols))
    val ret = mb.zeros(numRows, numCols)
    mm.addInPlace(false, this, that, ret)
    ret
  }

  def -[Repr1[C, D] >: Repr[C, D], Repr2[C, D] <: Matrix[C, D, Repr2], Repr3[C, D] <: Matrix[C, D, Repr3]]
        (that: Repr2[A, Raw])
        (implicit mm: MatrixAdder[A, Raw, Repr1, Repr2, Repr3],
                  mb: MatrixBuilder[A, Raw, Repr3]): Repr3[A, Raw] = {
    require(numRows == that.numRows && numCols == that.numCols,
            "Can't subtract matrices of dimensions (%d, %d) and (%d, %d)".format(numRows, numCols, that.numRows, that.numCols))
    val ret = mb.zeros(numRows, numCols)
    mm.addInPlace(true, this, that, ret)
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
trait MatrixAdder[A, Raw, Repr1[_, _], Repr2[_, _], Repr3[_, _]] {
  def addInPlace(sub: Boolean, m1: Repr1[A, Raw], m2: Repr2[A, Raw], ret: Repr3[A, Raw])
}

object MatrixMultiplier extends DenseMultipliers
trait MatrixMultiplier[A, Raw, Repr1[_, _], Repr2[_, _], Repr3[_, _]] {
  def gemm(alpha: A, beta: A, m1: Repr1[A, Raw], m2: Repr2[A, Raw], ret: Repr3[A, Raw])
}

object MatrixBuilder extends DenseBuilders
trait MatrixBuilder[A, Raw, Repr[_, _]] {
  def zeros(numRows: Int, numCols: Int): Repr[A, Raw]
}


object Test extends App {
  val m1 = MatrixBuilder.denseRow[Float, Float].zeros(4, 4)
  val m2 = MatrixBuilder.denseRow[Float, Float].zeros(4, 4)
  val m3: DenseRow[Float, Float] = m1+m2
  println(m3)
  val m4 = MatrixBuilder.denseRow[Double, Double].zeros(4, 4)
  m3.mapTo(_.toDouble, m4)

//  println(implicitly[MatrixAdder[Float, Float, DenseRow, DenseRow, Dense]])
}




/*
// why compiles?
trait Matrix[-Repr <: Matrix[Repr]]
trait Dense extends Matrix[Dense]
trait DenseRow extends Matrix[DenseRow] with Dense
*/