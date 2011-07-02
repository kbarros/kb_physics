package kip.math.fft

import kip.fftw3.{FFTW3Library => FFTW}
import FFTW.{INSTANCE => fft}

import com.sun.jna._
import java.nio.IntBuffer

object FFTReal {
  
  def convolve(dim: Array[Int], a: Array[Double], b: Array[Double], dst: Array[Double]) {
    val n = dim.product
    assert(a.size == n && b.size == n && dst.size == n)
    val fft = new FFTReal(dim)
    val ap = fft.allocFourierArray()
    val bp = fft.allocFourierArray()
    fft.forwardTransform(a, ap)
    fft.forwardTransform(b, bp)
    //fft.conjugateFourierArray(bp, bp)
    fft.multiplyFourierArrays(ap, bp, ap)
    fft.backwardTransform(ap, dst)
  }
  
  def test() {
    /*
    import kip.math.fft._
    */
    val dim = Array(3)
    val a = Array[Double](1,2,3)
    val b = Array[Double](1,2,1)
    val dst = new Array[Double](3)
    FFTReal.convolve(dim, a, b, dst)
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
  
  val in = fft.fftw_malloc(new NativeLong(inBytes))
  val out = fft.fftw_malloc(new NativeLong(outBytes))
  val inbuf = in.getByteBuffer(0, inBytes).asDoubleBuffer()
  val outbuf = out.getByteBuffer(0, outBytes).asDoubleBuffer()
  
  val planForward  = fft.fftw_plan_dft_r2c(dim.size, IntBuffer.wrap(dim), inbuf, outbuf, flags)
  val planBackward = fft.fftw_plan_dft_c2r(dim.size, IntBuffer.wrap(dim), outbuf, inbuf, flags)

  def forwardTransform(src: Array[Double], dst: Array[Double]) {
    assert(src.size == n)
    assert(dst.size == nrecip)
    
    inbuf.clear()
    inbuf.put(src)
    fft.fftw_execute(planForward)
    outbuf.rewind()
    outbuf.get(dst)
    
    // continuum normalization: f(k) = \int dx^d f(x) e^(i k x)
    val scale = len.product / dim.product
    for (i <- dst.indices) dst(i) *= scale
  }
  
  def backwardTransform(src: Array[Double], dst: Array[Double]) {
    assert(src.size == nrecip)
    assert(dst.size == n)
    
    outbuf.clear()
    outbuf.put(src)
    fft.fftw_execute(planBackward)
    inbuf.rewind()
    inbuf.get(dst)

    // continuum normalization: f(x) = (2 Pi)^(-d) \int dk^d f(k) e^(- i k x)
    val scale = 1 / len.product
    for (i <- dst.indices) dst(i) *= scale
  }
  
  def allocFourierArray(): Array[Double] = {
    new Array[Double](nrecip)
  }
  
  def multiplyFourierArrays(src1: Array[Double], src2: Array[Double], dst: Array[Double]) {
    assert(src1.size == nrecip)
    assert(src2.size == nrecip)
    assert(dst.size == nrecip)
    for (i <- 0 until src1.size/2) {
      // src and dst arrays might be aliased; create temporary variables
      val re = src1(2*i+0)*src2(2*i+0) - src1(2*i+1)*src2(2*i+1)
      val im = src1(2*i+0)*src2(2*i+1) + src1(2*i+1)*src2(2*i+0)
      dst(2*i+0) = re
      dst(2*i+1) = im
    }
  }

  def conjugateFourierArray(src: Array[Double], dst: Array[Double]) {
    assert(src.size == nrecip)
    assert(dst.size == nrecip)
    for (i <- 0 until src.size/2) {
      dst(2*i+0) = src(2*i+0)
      dst(2*i+1) = -src(2*i+1)
    }
  }
  
  // for each indexed complex number in fourier array, return corresponding vector k
  // where component k_r = n (2 pi / L_r) for integer n in range [-N/2, +N/2)
  def fourierVector(i: Int): Array[Double] = {
    assert (0 <= i && i < nrecip)
    
    val k = new Array[Double](rank)
    var ip = i
    for (r <- rank-1 to 0 by -1) {
      val d = if (r == rank-1) (dim(r)/2+1) else dim(r) // fftw compresses last index
      k(r) = ip % d
      if (k(r) >= d/2)
        k(r) -= d
      val dk = 2*math.Pi/len(r)
      k(r) *= dk
      ip /= d
    }
    
    k
  }
  
  def destroy {
    fft.fftw_destroy_plan(planForward)
    fft.fftw_destroy_plan(planBackward)
    fft.fftw_free(in)
    fft.fftw_free(out)
  }

}
