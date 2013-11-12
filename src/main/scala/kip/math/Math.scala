package kip.math


object Math {
  def sqr(x: Double) = x*x
  def cube(x: Double) = x*x*x
  
  def mod(x: Double, y: Double) = ((x % y) + y) % y
  def mod(x: Int, y: Int) = ((x % y) + y) % y
  
  def fftReal(in: Array[Double]): Array[Complex] = {
    val tran = new fft.FFTReal(dim=Array(in.size))
    val recipArray = tran.allocFourierArray()
    tran.forwardTransform(src=in, dst=recipArray)
    recipArray.grouped(2).map{a => Complex(a(0), a(1))}.toArray
  }
}
