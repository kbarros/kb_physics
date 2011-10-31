package kip.util


object Statistics {

  class OnlineVariance {
    private var _n = 0.0
    private var _mean = 0.0
    private var m2 = 0.0 // second central moment, sum(i <- 1 to n) (x_i - mean_n)^2
    
    def n = _n
    def mean = _mean
    def variance = m2 / n
    def sample_variance = m2 / (n-1) // sample variance; an unbiased estimator of true variance
    def stddev = math.sqrt(sample_variance) // sample standard deviation
    
    def accum(x: Double) = {
      _n += 1.0
      val delta = x - mean
      _mean += delta/n
      m2 += delta*(x - mean)
    }
  }

  class OnlineMoments {
    private var _n = 0.0
    private var _mean = 0.0
    private var m2, m3, m4 = 0.0 // higher order central moments
    
    def n = _n
    def mean = _mean
    def variance = m2 / n
    def sample_variance = m2 / (n-1) // sample variance; an unbiased estimator of true variance
    def stddev = math.sqrt(sample_variance) // sample standard deviation
    def skewness = (m3/n) / math.pow(m2/n, 3./2.)
    def kurtosis = (n*m4) / (m2*m2) - 3
    
    def accum(x: Double) = {
      _n += 1
      val delta = x - mean
      val delta_n = delta / n
      val delta_n2 = delta_n * delta_n
      val term1 = delta * delta_n * (n-1)
      _mean += delta_n
      m4 += term1 * delta_n2 * (n*n - 3*n + 3) + 6 * delta_n2 * m2 - 4 * delta_n * m3
      m3 += term1 * delta_n * (n - 2) - 3 * delta_n * m2
      m2 += term1
    }
  }
  
  
  def mean(vs: Traversable[Double]) = {
    if (vs.size == 0) Double.NaN else vs.sum / vs.size
  }
  
  // Approximate error in the mean
  def mean_err(vs: Traversable[Double]) = {
    stddev(vs) / math.sqrt(vs.size)
  }
  
  def stddev(vs: Traversable[Double]) = {
    val o = new OnlineVariance
    vs.foreach(o accum _)
    o.stddev
  }
  
  def skewness(vs: Traversable[Double]) = {
    val o = new OnlineMoments
    vs.foreach(o accum _)
    o.skewness
  }

  def kurtosis(vs: Traversable[Double]) = {
    val o = new OnlineMoments
    vs.foreach(o accum _)
    o.kurtosis
  }
}
