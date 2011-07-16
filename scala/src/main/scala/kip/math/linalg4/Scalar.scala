package kip.math.linalg4

import kip.math.Complex

import java.nio.{DoubleBuffer, FloatBuffer}

object Scalar {
  trait RealDbl extends Scalar { type A = Double; type Raw = Double; type Buf = DoubleBuffer }
  trait RealFlt extends Scalar { type A = Float; type Raw = Float; type Buf = FloatBuffer }
  trait ComplexDbl extends Scalar { type A = Complex; type Raw = Double; type Buf = DoubleBuffer }
  trait ComplexFlt extends Scalar { type A = Complex; type Raw = Float; type Buf = FloatBuffer }
}

trait Scalar {
  type A
  type Raw
  type Buf
}


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
  
  object Builder {
    implicit val RealFlt = new Builder[Scalar.RealFlt] {
      def build(size: Int) = new RealFlt(size)
    }
    implicit val RealDbl = new Builder[Scalar.RealDbl] {
      def build(size: Int) = new RealDbl(size)
    }
    implicit val ComplexFlt = new Builder[Scalar.ComplexFlt] {
      def build(size: Int) = new ComplexFlt(size)
    }
    implicit val ComplexDbl = new Builder[Scalar.ComplexDbl] {
      def build(size: Int) = new ComplexDbl(size)
    }
  }
  
  trait Builder[S <: Scalar] {
    def build(size: Int): ScalarData[S]
  }  


  class RealFlt(size: Int) extends ScalarData[Scalar.RealFlt] {
    val scalar = ScalarOps.RealFlt

    val raw = new Array[Float](size)
    val buffer = java.nio.FloatBuffer.wrap(raw)
    
    def rawApply(i: Int): Float = raw(i)
    def rawUpdate(i: Int, x: Float) { raw(i) = x }
    
    def apply(i: Int): Float = raw(i)
    def update(i: Int, x: Float) = raw(i) = x
    override def madd(i0: Int, a1: ScalarData[Scalar.RealFlt], i1: Int, a2: ScalarData[Scalar.RealFlt], i2: Int) {
      raw(i0) += a1.rawApply(i1)*a2.rawApply(i2)
    }
  }

  class RealDbl(size: Int) extends ScalarData[Scalar.RealDbl] {
    val scalar = ScalarOps.RealDbl
    
    val raw = new Array[Double](size)
    val buffer = java.nio.DoubleBuffer.wrap(raw)
    
    def rawApply(i: Int): Double = raw(i)
    def rawUpdate(i: Int, x: Double) { raw(i) = x }
    
    def apply(i: Int): Double = raw(i)
    def update(i: Int, x: Double) = raw(i) = x
    override def madd(i0: Int, a1: ScalarData[Scalar.RealDbl], i1: Int, a2: ScalarData[Scalar.RealDbl], i2: Int) {
      raw(i0) += a1.rawApply(i1)*a2.rawApply(i2)
    }
  }

  class ComplexFlt(size: Int) extends ScalarData[Scalar.ComplexFlt] {
    val scalar = ScalarOps.ComplexDbl
    
    val raw = new Array[Float](2*size)
    val buffer = java.nio.FloatBuffer.wrap(raw)

    def rawApply(i: Int): Float = raw(i)
    def rawUpdate(i: Int, x: Float) { raw(i) = x }
    
    def apply(i: Int) = Complex(raw(2*i+0), raw(2*i+1))
    def update(i: Int, x: Complex) {
      raw(2*i+0) = x.re.toFloat
      raw(2*i+1) = x.im.toFloat
    }
    override def madd(i0: Int, a1: ScalarData[Scalar.ComplexFlt], i1: Int, a2: ScalarData[Scalar.ComplexFlt], i2: Int) {
      val x1_re = a1.rawApply(2*i1+0)
      val x1_im = a1.rawApply(2*i1+1)
      val x2_re = a2.rawApply(2*i2+0)
      val x2_im = a2.rawApply(2*i2+1)
      raw(2*i0+0) += x1_re*x2_re - x1_im*x2_im
      raw(2*i0+1) += x1_re*x2_im + x1_im*x2_re
    }
  }

  class ComplexDbl(size: Int) extends ScalarData[Scalar.ComplexDbl] {
    val scalar = ScalarOps.ComplexDbl
    
    val raw = new Array[Double](2*size)
    val buffer = java.nio.DoubleBuffer.wrap(raw)

    def rawApply(i: Int): Double = raw(i)
    def rawUpdate(i: Int, x: Double) { raw(i) = x }
    
    def apply(i: Int) = Complex(raw(2*i+0), raw(2*i+1))
    def update(i: Int, x: Complex) {
      raw(2*i+0) = x.re.toDouble
      raw(2*i+1) = x.im.toDouble
    }
    override def madd(i0: Int, a1: ScalarData[Scalar.ComplexDbl], i1: Int, a2: ScalarData[Scalar.ComplexDbl], i2: Int) {
      val x1_re = a1.rawApply(2*i1+0)
      val x1_im = a1.rawApply(2*i1+1)
      val x2_re = a2.rawApply(2*i2+0)
      val x2_im = a2.rawApply(2*i2+1)
      raw(2*i0+0) += x1_re*x2_re - x1_im*x2_im
      raw(2*i0+1) += x1_re*x2_im + x1_im*x2_re
    }
  }
}


trait ScalarData[S <: Scalar] {
  def scalar: ScalarOps[S#A]
  def buffer: S#Buf
  def rawApply(i: Int): S#Raw
  def rawUpdate(i: Int, x: S#Raw)
  def apply(i: Int): S#A
  def update(i: Int, x: S#A)
  def madd(i0: Int, a1: ScalarData[S], i1: Int, a2: ScalarData[S], i2: Int) {
    this(i0) = scalar.add(this(i0), scalar.mul(a1(i1), a2(i2)))
  }
  def dispose() = ()
}
