package kip.math
package linalg


object ScalarData {
  implicit object RealDbl extends RealDblData
  implicit object RealFlt extends RealFltData
  implicit object ComplexDbl extends ComplexDblData
  implicit object ComplexFlt extends ComplexFltData
}

trait ScalarData[@specialized(Float, Double) A, B] extends Scalar[A] {
  type Buf
  def alloc(size: Int): B
  def buffer(a: B): Buf
  def copy(src: B, srcPos: Int, dst: B, dstPos: Int, length: Int): Unit
  def update(a: B, i: Int, x: A): Unit
  def apply(a: B, i: Int): A
  def size(a: B): Int
  def dispose(a: B): Unit
  
  def clone(a: B): B = {
    val ret = alloc(size(a))
    copy(a, 0, ret, 0, size(a))
    ret
  }
    
  def madd(a0: B, i0: Int, a1: B, i1: Int, a2: B, i2: Int) {
    val x0 = apply(a0, i0)
    val x1 = apply(a1, i1)
    val x2 = apply(a2, i2)
    update(a0, i0, add(x0, mul(x1, x2)))
  }
  
  def madd(a0: B, i0: Int, a1: B, i1: Int, x2: A) {
    val x0 = apply(a0, i0)
    val x1 = apply(a1, i1)
    update(a0, i0, add(x0, mul(x1, x2)))
  }
}


// ------------------------------------------------------------
// Element types of packaged buffer

trait ScalarBufferedDbl {
  type Buf = Array[Double] // java.nio.DoubleBuffer
}

trait ScalarBufferedFlt {
  type Buf = Array[Float] // java.nio.FloatBuffer
}


// ------------------------------------------------------------
// Real

trait RealDblData extends ScalarData[Double, Array[Double]] with ScalarBufferedDbl with RealDblTC {
  def alloc(size: Int) = new Array[Double](size)
  def buffer(a: Array[Double]) = a
  def copy(src: Array[Double], srcPos: Int, dst: Array[Double], dstPos: Int, length: Int) {
    System.arraycopy(src, srcPos, dst, dstPos, length)
  }
  def update(a: Array[Double], i: Int, x: Double) { a(i) = x }
  def apply(a: Array[Double], i: Int): Double = a(i)
  def size(a: Array[Double]) = a.size
  def dispose(a: Array[Double]) { }
}

trait RealFltData extends ScalarData[Float, Array[Float]] with ScalarBufferedFlt with RealFltTC {
  def alloc(size: Int) = new Array[Float](size)
  def buffer(a: Array[Float]) = a
  def copy(src: Array[Float], srcPos: Int, dst: Array[Float], dstPos: Int, length: Int) {
    System.arraycopy(src, srcPos, dst, dstPos, length)
  }
  def update(a: Array[Float], i: Int, x: Float) { a(i) = x }
  def apply(a: Array[Float], i: Int): Float = a(i)
  def size(a: Array[Float]) = a.size
  def dispose(a: Array[Float]) { }
}


// ------------------------------------------------------------
// Complex

trait ComplexDblData extends ScalarData[Complex, Array[Double]] with ScalarBufferedDbl with ComplexTC {
  def alloc(size: Int) = new Array[Double](2*size)
  def buffer(a: Array[Double]) = a
  def copy(src: Array[Double], srcPos: Int, dst: Array[Double], dstPos: Int, length: Int) {
    System.arraycopy(src, 2*srcPos, dst, 2*dstPos, 2*length)
  }
  def update(a: Array[Double], i: Int, x: Complex) {
    a(2*i+0) = x.re
    a(2*i+1) = x.im
  }
  def apply(a: Array[Double], i: Int): Complex = Complex(a(2*i+0), a(2*i+1))
  def size(a: Array[Double]) = a.size/2
  def dispose(a: Array[Double]) { }
  override def madd(a0: Array[Double], i0: Int, a1: Array[Double], i1: Int, a2: Array[Double], i2: Int) {
    val a1_re = a1(2*i1+0)
    val a1_im = a1(2*i1+1)
    val a2_re = a2(2*i2+0)
    val a2_im = a2(2*i2+1)
    a0(2*i0+0) += a1_re*a2_re - a1_im*a2_im
    a0(2*i0+1) += a1_re*a2_im + a1_im*a2_re
  }
  override def madd(a0: Array[Double], i0: Int, a1: Array[Double], i1: Int, x2: Complex) {
    val a1_re = a1(2*i1+0)
    val a1_im = a1(2*i1+1)
    a0(2*i0+0) += a1_re*x2.re - a1_im*x2.im
    a0(2*i0+1) += a1_re*x2.im + a1_im*x2.re
  }
}

trait ComplexFltData extends ScalarData[Complex, Array[Float]] with ScalarBufferedFlt with ComplexTC {
  def alloc(size: Int) = new Array[Float](2*size)
  def buffer(a: Array[Float]) = a
  def copy(src: Array[Float], srcPos: Int, dst: Array[Float], dstPos: Int, length: Int) {
    System.arraycopy(src, 2*srcPos, dst, 2*dstPos, 2*length)
  }
  def update(a: Array[Float], i: Int, x: Complex) {
    a(2*i+0) = x.re.toFloat
    a(2*i+1) = x.im.toFloat
  }
  def apply(a: Array[Float], i: Int): Complex = Complex(a(2*i+0), a(2*i+1))
  def size(a: Array[Float]) = a.size/2
  def dispose(a: Array[Float]) { }
  override def madd(a0: Array[Float], i0: Int, a1: Array[Float], i1: Int, a2: Array[Float], i2: Int) {
    val a1_re = a1(2*i1+0)
    val a1_im = a1(2*i1+1)
    val a2_re = a2(2*i2+0)
    val a2_im = a2(2*i2+1)
    a0(2*i0+0) += a1_re*a2_re - a1_im*a2_im
    a0(2*i0+1) += a1_re*a2_im + a1_im*a2_re
  }
  override def madd(a0: Array[Float], i0: Int, a1: Array[Float], i1: Int, x2: Complex) {
    val a1_re = a1(2*i1+0)
    val a1_im = a1(2*i1+1)
    a0(2*i0+0) += (a1_re*x2.re - a1_im*x2.im).toFloat
    a0(2*i0+1) += (a1_re*x2.im + a1_im*x2.re).toFloat
  }
}


