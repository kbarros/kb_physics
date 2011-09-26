package kip.math


object Math {
  def sqr(x: Double) = x*x
  def hypot(x: Double, y: Double, z: Double) = math.sqrt(x*x + y*y + z*z)
  def avg(x: Double, y: Double) = 0.5*(x+y)
  
  def fftReal(in: Array[Double]): Array[Complex] = {
    val tran = new fft.FFTReal(dim=Array(in.size))
    val recipArray = tran.allocFourierArray()
    tran.forwardTransform(src=in, dst=recipArray)
    recipArray.grouped(2) map { a => Complex(a(0), a(1)) } toArray
  }
}
