package kip.math.fft

import com.sun.jna._
import java.nio.IntBuffer
import kip.fftw3.{FFTW3Library => FFTW}
import FFTW.{INSTANCE => fftw}


object FFTComplex {

}


// Note that arrays are packed in row major order, so the last index is the fastest varying.
// Thus, if indices are computed as (i = Lx*y + x) then one should use dim(Ly, Lx)
class FFTComplex(dim: Array[Int], lenOption: Option[Array[Double]] = None, flags: Int = FFTW.FFTW_ESTIMATE) {

  val len = lenOption.getOrElse(dim.map(_.toDouble))
  val rank = dim.size
  
  // number of doubles in real/reciprocal-space array
  val n = 2*dim.product
  
  val sizeofDouble = 8;
  val inBytes  = sizeofDouble*n
  val outBytes = sizeofDouble*n
  
  val in = fftw.fftw_malloc(new NativeLong(inBytes))
  val out = fftw.fftw_malloc(new NativeLong(outBytes))
  val inbuf = in.getByteBuffer(0, inBytes).asDoubleBuffer()
  val outbuf = out.getByteBuffer(0, outBytes).asDoubleBuffer()
  
  val planForward  = fftw.fftw_plan_dft(rank, IntBuffer.wrap(dim), inbuf, outbuf, FFTW.FFTW_FORWARD, flags)
  val planBackward = fftw.fftw_plan_dft(rank, IntBuffer.wrap(dim), outbuf, inbuf, FFTW.FFTW_BACKWARD, flags)

  def forwardTransform(src: Array[Double], dst: Array[Double]) {
    require(src.size == n)
    require(dst.size == n)
    
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
    require(src.size == n)
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
    new Array[Double](n)
  }
  
  def tabulateFourierArray(dst: Array[Double])(f: Array[Double] => (Double, Double)) {
    require(dst.size == n)
    for (i <- 0 until dst.size/2) {
      val k = fourierVector(i)
      val (re, im) = f(k)
      dst(2*i+0) = re
      dst(2*i+1) = im
    }
  }
    
  // Returns the list of all fourier vectors
  def fourierVectors: Array[Array[Double]] = {
    Array.tabulate(n/2) { fourierVector(_) }
  }
  
  // for each indexed complex number in fourier array, return corresponding vector k
  // where component k(r) = n (2 pi / L_r) for integer n in range [-N/2, +N/2)
  def fourierVector(i: Int): Array[Double] = {
    require(0 <= i && i < n/2)
    val k = new Array[Double](rank)
    var ip = i
    for (r <- rank-1 to 0 by -1) {
      k(r) = ip % dim(r)
      if (k(r) >= dim(r)/2)
        k(r) -= dim(r)
      val dk = 2*math.Pi/len(r)
      k(r) *= dk
      ip /= dim(r)
    }
    k
  }
  
  def destroy {
    fftw.fftw_destroy_plan(planForward)
    fftw.fftw_destroy_plan(planBackward)
    fftw.fftw_free(in)
    fftw.fftw_free(out)
  }

}