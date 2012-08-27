package kip.math

package object fft {
  def multiplyComplexArrays(src1: Array[Double], src2: Array[Double], dst: Array[Double]) {
    require(src1.size == src2.size && src2.size == dst.size)
    for (i <- 0 until src1.size/2) {
      // src and dst arrays might be aliased; create temporary variables
      val re = src1(2*i+0)*src2(2*i+0) - src1(2*i+1)*src2(2*i+1)
      val im = src1(2*i+0)*src2(2*i+1) + src1(2*i+1)*src2(2*i+0)
      dst(2*i+0) = re
      dst(2*i+1) = im
    }
  }

  def conjugateComplexArray(src: Array[Double], dst: Array[Double]) {
    require(src.size == dst.size)
    for (i <- 0 until src.size/2) {
      dst(2*i+0) = src(2*i+0)
      dst(2*i+1) = -src(2*i+1)
    }
  }

}