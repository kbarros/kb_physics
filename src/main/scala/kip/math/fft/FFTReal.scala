package kip.math.fft

import com.sun.jna._
import java.nio.IntBuffer
import kip.fftw3.{FFTW3Library => FFTW}
import FFTW.{INSTANCE => fftw}

object FFTReal {
  
  def test() {
    /*
    import kip.math.fft._
    */
    val dim = Array(3, 3)
    val a = Array[Double](0,0,1, 0,0,0, 0,0,0)
    val b = Array[Double](1,2,3, 4,5,6, 7,8,9)
    val dst = new Array[Double](dim.product)
    val fft = new FFTReal(dim)
    fft.convolve(a, b, dst)
    dst.foreach(println _)
  }

  def test2() {
    val dim = Array(4, 4)
    val a = Array[Double](0,1,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0)
    val fft = new FFTReal(dim)
    val ap = fft.allocFourierArray()
    fft.forwardTransform(a, ap)
    
    val bp = fft.allocFourierArray()
    for (i <- 0 until bp.size/2) {
      val k = fft.fourierVector(i)
      bp(2*i+0) = math.cos(-k(1))
      bp(2*i+1) = math.sin(-k(1))
    }
    
    for (i <- ap.indices) {
      println(ap(i) - bp(i))
    }
  }

  def test3() {
    val dim = Array(4, 4)
    val a = Array[Double](0,0,0,0, 1,0,0,0, 0,0,0,0, 0,0,0,0)
    val fft = new FFTReal(dim)
    val ap = fft.allocFourierArray()
    fft.forwardTransform(a, ap)
    
    val bp = fft.allocFourierArray()
    for (i <- 0 until bp.size/2) {
      val k = fft.fourierVector(i)
      bp(2*i+0) = math.cos(-k(0))
      bp(2*i+1) = math.sin(-k(0))
    }
    
    for (i <- ap.indices) {
      println(ap(i) - bp(i))
    }
  }
}


class FFTReal(dim: Array[Int], lenOption: Option[Array[Double]] = None, flags: Int = FFTW.FFTW_ESTIMATE) {

  val len = lenOption.getOrElse(dim.map(_.toDouble))
  val rank = dim.size
  
  // number of doubles in real-space array
  val n = dim.product
  // number of doubles in reciprocal space array
  val nrecip = 2*(dim.slice(0,rank-1).product*(dim(rank-1)/2+1)) // fftw compresses last index
  
  val sizeofDouble = 8;
  val inBytes  = sizeofDouble*n
  val outBytes = sizeofDouble*nrecip
  
  val in = fftw.fftw_malloc(new NativeLong(inBytes))
  val out = fftw.fftw_malloc(new NativeLong(outBytes))
  val inbuf = in.getByteBuffer(0, inBytes).asDoubleBuffer()
  val outbuf = out.getByteBuffer(0, outBytes).asDoubleBuffer()
  
  val planForward  = fftw.fftw_plan_dft_r2c(dim.size, IntBuffer.wrap(dim), inbuf, outbuf, flags)
  val planBackward = fftw.fftw_plan_dft_c2r(dim.size, IntBuffer.wrap(dim), outbuf, inbuf, flags)

  def forwardTransform(src: Array[Double], dst: Array[Double]) {
    require(src.size == n)
    require(dst.size == nrecip)
    
    inbuf.clear()
    inbuf.put(src)
    fftw.fftw_execute(planForward)
    outbuf.rewind()
    outbuf.get(dst)
    
    // continuum normalization: f(k) = \int dx^d f(x) e^(i k x)
    val scale = len.product / dim.product
    for (i <- dst.indices) dst(i) *= scale
  }
  
  def backwardTransform(src: Array[Double], dst: Array[Double]) {
    require(src.size == nrecip)
    require(dst.size == n)
    
    outbuf.clear()
    outbuf.put(src)
    fftw.fftw_execute(planBackward)
    inbuf.rewind()
    inbuf.get(dst)

    // continuum normalization: f(x) = (2 Pi)^(-d) \int dk^d f(k) e^(- i k x)
    val scale = 1 / len.product
    for (i <- dst.indices) dst(i) *= scale
  }
  
  def allocFourierArray(): Array[Double] = {
    new Array[Double](nrecip)
  }
  
  // TODO: generalize to arbitrary rank
  def uncompressFourierArray(f: Array[Double]): Array[Double] = {
    val ret = new Array[Double](2*dim.product)
    require(rank == 2)
    for (k0 <- 0 until dim(0); k1 <- 0 until dim(1)) {
      val i = k1 + dim(1)*(k0)
      val j = arrayIndex(Array(k0, k1))
      if (j >= 0) {
        ret(2*i+0) = f(2*j+0)
        ret(2*i+1) = f(2*j+1)
      }
      else {
        ret(2*i+0) = + f(-2*j+0)
        ret(2*i+1) = - f(-2*j+1) // complex conjugate
      }
    }
    ret
  }

  lazy val allFourierVectors = Array.tabulate(nrecip/2)(fourierVector _)
  
  // for each indexed complex number in fourier array, return corresponding vector k
  // where component k(r) = n (2 pi / L_r) for integer n in range [-N/2, +N/2)
  def fourierVector(i: Int): Array[Double] = {
    val k = fourierIndices(i)
    Array.tabulate[Double](rank) (r => 2*math.Pi*k(r) / len(r))
  }
  
  def fourierIndices(i: Int): Array[Int] = {
    require(0 <= i && i < nrecip/2)
    val k = new Array[Int](rank)
    var ip = i
    for (r <- rank-1 to 0 by -1) {
      val d = if (r == rank-1) (dim(r)/2+1) else dim(r) // fftw compresses last index
      k(r) = ip % d
      if (k(r) >= dim(r)/2)
        k(r) -= dim(r)
      ip /= d
    }
    k
  }
  
  // Returns the array index corresponding to Fourier vector indices
  //
  // Note that arrays are packed in row major order, so the last index is the fastest varying.
  // For example, in three dimensions: 
  //  system dimensions   = (Lx, Ly, Lz)
  //  fourier k-indices   = (kx, ky, kz)
  //  fourier array index = kz + Lz * (ky + ly * (kz))
  //
  // Not all Fourier vectors k map to valid array indices. In such cases, we encode ki < 0
  // such that: f(ki) == f(-ki)^\ast
  def arrayIndex(k: Array[Int]): Int = {
    import kip.math.Math.mod
    
    // no index exists for some fourier indices. these are available by complex conjugation
    // of k -> -k index
    val invert = if (mod(k(rank-1), dim(rank-1)) > dim(rank-1)/2) -1 else 1
    
    val kp = Array.tabulate[Int](rank) { r => mod(invert * k(r), dim(r)) }
    
    var ip = 0
    for (r <- 0 until rank) {
      val d = if (r == rank-1) (dim(r)/2+1) else dim(r)
      ip *= d
      ip += kp(r)
    }
    
    invert * ip
  }

  def destroy {
    fftw.fftw_destroy_plan(planForward)
    fftw.fftw_destroy_plan(planBackward)
    fftw.fftw_free(in)
    fftw.fftw_free(out)
  }

  
  // returns dst(i) = \sum_j a(i) b(j-i)
  def convolve(a: Array[Double], b: Array[Double], dst: Array[Double]) {
    require(a.size == n && b.size == n && dst.size == n)
    val ap = allocFourierArray()
    val bp = allocFourierArray()
    forwardTransform(a, ap)
    forwardTransform(b, bp)
    multiplyComplexArrays(ap, bp, ap)
    backwardTransform(ap, dst)
  }

  def convolveWithRecip(a: Array[Double], bp: Array[Double], dst: Array[Double]) {
    require(a.size == n && bp.size == nrecip && dst.size == n)
    val ap = allocFourierArray()
    forwardTransform(a, ap)
    multiplyComplexArrays(ap, bp, ap)
    backwardTransform(ap, dst)
  }
}
