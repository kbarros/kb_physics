package kip.math.linalg4

import kip.math.Complex
import java.nio.{DoubleBuffer, FloatBuffer}
import java.nio.Buffer

object Scalar {
  trait RealTyp    extends Scalar
  trait ComplexTyp extends Scalar
  
  trait RealDbl    extends RealTyp    { type A = Double;  type Raw = Double; type Buf = DoubleBuffer }
  trait RealFlt    extends RealTyp    { type A = Float;   type Raw = Float;  type Buf = FloatBuffer }
  trait ComplexDbl extends ComplexTyp { type A = Complex; type Raw = Double; type Buf = DoubleBuffer }
  trait ComplexFlt extends ComplexTyp { type A = Complex; type Raw = Float;  type Buf = FloatBuffer }
}

trait Scalar {
  type A
  type Raw
  type Buf <: Buffer
}


object ScalarOps {
  implicit object RealFlt extends RealFlt
  implicit object RealDbl extends RealDbl
  implicit object ComplexFlt extends ComplexFlt
  implicit object ComplexDbl extends ComplexDbl

  trait RealFlt extends ScalarOps[Scalar.RealFlt] {
    def add(a: Float, b: Float): Float = a + b
    def sub(a: Float, b: Float): Float = a - b
    def mul(a: Float, b: Float): Float = a * b
    def div(a: Float, b: Float): Float = a / b
    def neg(a: Float): Float = -a
    def conj(a: Float): Float = a
    def zero: Float = 0.0f
    def one: Float = 1.0f
    
    def components = 1
    def read(data: Data, i: Int): Float = data(i)
    def write(data: Data, i: Int, x: Float): Unit = data(i) = x
    override def madd(data0: Data, i0: Int, data1: Data, i1: Int, x2: Float) {
      data0(i0) += data1(i1) * x2
    }
    override def madd(data0: Data, i0: Int, data1: Data, i1: Int, data2: Data, i2: Int) {
      data0(i0) += data1(i1) * data2(i2)
    }
  }

  trait RealDbl extends ScalarOps[Scalar.RealDbl] {
    def add(a: Double, b: Double): Double = a + b
    def sub(a: Double, b: Double): Double = a - b
    def mul(a: Double, b: Double): Double = a * b
    def div(a: Double, b: Double): Double = a / b
    def neg(a: Double): Double = -a
    def conj(a: Double): Double = a
    def zero: Double = 0.0
    def one: Double = 1.0
    
    def components = 1
    def read(data: Data, i: Int): Double = data(i)
    def write(data: Data, i: Int, x: Double): Unit = data(i) = x
    override def madd(data0: Data, i0: Int, data1: Data, i1: Int, x2: Double) {
      data0(i0) += data1(i1) * x2
    }
    override def madd(data0: Data, i0: Int, data1: Data, i1: Int, data2: Data, i2: Int) {
      data0(i0) += data1(i1) * data2(i2)
    }
  }

  trait ComplexFlt extends ScalarOps[Scalar.ComplexFlt] {
    def add(a: Complex, b: Complex): Complex = a + b
    def sub(a: Complex, b: Complex): Complex = a - b
    def mul(a: Complex, b: Complex): Complex = a * b
    def div(a: Complex, b: Complex): Complex = a / b
    def neg(a: Complex): Complex = -a
    def conj(a: Complex): Complex = a.conj
    def zero: Complex = 0.0f
    def one: Complex = 1.0f
    
    def components = 2
    def read(data: Data, i: Int) = Complex(data(2*i+0), data(2*i+1))
    def write(data: Data, i: Int, x: Complex) {
      data(2*i+0) = x.re.toFloat
      data(2*i+1) = x.im.toFloat
    }
    override def madd(data0: Data, i0: Int, data1: Data, i1: Int, x2: Complex) {
      val x1_re = data1(2*i1+0)
      val x1_im = data1(2*i1+1)
      data0(2*i0+0) += x1_re*x2.re.toFloat - x1_im*x2.im.toFloat
      data0(2*i0+1) += x1_re*x2.im.toFloat + x1_im*x2.re.toFloat
    }
    override def madd(data0: Data, i0: Int, data1: Data, i1: Int, data2: Data, i2: Int) {
      val x1_re = data1(2*i1+0)
      val x1_im = data1(2*i1+1)
      val x2_re = data2(2*i2+0)
      val x2_im = data2(2*i2+1)
      data0(2*i0+0) += x1_re*x2_re - x1_im*x2_im
      data0(2*i0+1) += x1_re*x2_im + x1_im*x2_re
    }
  }

  trait ComplexDbl extends ScalarOps[Scalar.ComplexDbl] {
    def add(a: Complex, b: Complex): Complex = a + b
    def sub(a: Complex, b: Complex): Complex = a - b
    def mul(a: Complex, b: Complex): Complex = a * b
    def div(a: Complex, b: Complex): Complex = a / b
    def neg(a: Complex): Complex = -a
    def conj(a: Complex): Complex = a.conj
    def zero: Complex = 0.0f
    def one: Complex = 1.0f
    
    def components = 2
    def read(data: Data, i: Int) = Complex(data(2*i+0), data(2*i+1))
    def write(data: Data, i: Int, x: Complex) {
      data(2*i+0) = x.re
      data(2*i+1) = x.im
    }
    override def madd(data0: Data, i0: Int, data1: Data, i1: Int, x2: Complex) {
      val x1_re = data1(2*i1+0)
      val x1_im = data1(2*i1+1)
      data0(2*i0+0) += x1_re*x2.re - x1_im*x2.im
      data0(2*i0+1) += x1_re*x2.im + x1_im*x2.re
    }
    override def madd(data0: Data, i0: Int, data1: Data, i1: Int, data2: Data, i2: Int) {
      val x1_re = data1(2*i1+0)
      val x1_im = data1(2*i1+1)
      val x2_re = data2(2*i2+0)
      val x2_im = data2(2*i2+1)
      data0(2*i0+0) += x1_re*x2_re - x1_im*x2_im
      data0(2*i0+1) += x1_re*x2_im + x1_im*x2_re
    }
  }
}


trait GenScalarOps[@specialized(Float, Double) A, @specialized(Float, Double) Raw, Buf <: Buffer] {
  def add(a: A, b: A): A
  def sub(a: A, b: A): A
  def mul(a: A, b: A): A
  def div(a: A, b: A): A
  def neg(a: A): A
  def conj(a: A): A
  def zero: A
  def one: A
  
  type Data = RawData[Raw, Buf]
  def components: Int
  def read(data: Data, i: Int): A
  def write(data: Data, i: Int, x: A)
  def madd(data0: Data, i0: Int, data1: Data, i1: Int, x2: A) {
    write(data0, i0, add(read(data0, i0), mul(read(data1, i1), x2)))
  }
  def madd(data0: Data, i0: Int, data1: Data, i1: Int, data2: Data, i2: Int) {
    write(data0, i0, add(read(data0, i0), mul(read(data1, i1), read(data2, i2))))
  }
}
trait ScalarOps[S <: Scalar] extends GenScalarOps[S#A, S#Raw, S#Buf]
