package kip.util


class BlockAnalysis(data: Array[Double]) {
  case class Block(sigma: Double, sigma_err: Double, corr: Double, n: Int)

  val mean = data.sum / data.size
  
  val blocks = {
    var a = data
    var n = a.size
    var ret: List[Block] = Nil
    
    while (n > 16) {
      val c0 = c_t(a, 0)
      val c1 = c_t(a, 1)

      // see Flyvbjerb & Petersen, 1989
      val sigma2 = c0 / (n - 1.0)
      val sigma2_err = sigma2 * math.sqrt(2.0 / (n - 1.0))
      
      // see Propagation of Uncertainty, Wikipedia,
      // (x = a +- e) and (y = x^c) => (y = a^c +- c e (a^c / a))
      val sigma = math.sqrt(sigma2)
      val sigma_err = 0.5 * sigma2_err * sigma / sigma2

      val corr = c1 / c0
      ret = Block(sigma, sigma_err, corr, n) +: ret

      a = decimate(a, 2)
      n = a.size
    }
    
    // finest blocking is first in return list
    ret.reverse
  }
  
  val (error, error_error, isDecorrelated) = {
    blocks.find(_.corr < 0.2) match {
      case Some(iter) => (iter.sigma, iter.sigma_err, true)
      case None       => (blocks.last.sigma, Double.PositiveInfinity, false)
    }
  }
  
  private def decimate(a: Array[Double], factor: Int): Array[Double] = {
    val n = a.size
    val ret = new Array[Double](n / factor)
    for (i <- ret.indices) {
      ret(i) = 0
      for (j <- 0 until factor) {
        ret(i) += a(factor*i+j)/factor
      }
    }
    ret
  }
  
  private def c_t(a: Array[Double], t: Int) = {
    val n = a.size
    var ret = 0.0
    for (i <- 0 until (n-t)) {
      val del1 = a(i+0) - mean
      val del2 = a(i+t) - mean
      ret += del1*del2
    }
    ret / (n-t)
  }
}
