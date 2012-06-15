package kip.projects.quantum

import scala.util.Random

// solves the equation
//   d^2x/dt^2 = - dF/dx - gamma \dot x + sqrt(2 T gamma) eta
//
// via
//
//   dx/dt = p
//   dp/dt = - dF/dx - gamma p + sqrt(2 T gamma) eta
//
// using a Verlet like algorithm (Trotter expansion in position and momentum operators)
//

abstract class InertialLangevin(x: Array[R], gamma: R, T: R, dt: R, subIter: Int, rand: Random) {
  val N = x.size
  val p = new Array[R](N)
  val f = new Array[R](N)
  
  // fills force array, f, given position x
  def calcForce(x: Array[R], f: Array[R])
  
  // projects vector x onto base manifold
  def projectToBase(x: Array[R]) { }
  
  // projects vector t (e.g., momentum, force) onto tangent manifold at point x
  def projectToTangent(x: Array[R], t: Array[R]) { }
  
  def step() {
    val dt2 = dt / subIter
    
    // update momentum at fixed position
    calcForce(x, f)
    for (k <- 0 until subIter) {
      for (i <- 0 until N) {
        p(i) = p(i) + dt2 * (- f(i) - gamma * p(i)) + (math.sqrt(dt2 * 2 * T * gamma) * rand.nextGaussian()).toFloat
      }
      projectToTangent(x, p)
    }
    
    // update position at fixed momentum
    for (k <- 0 until subIter) {
      for (i <- 0 until N) {
        x(i) = x(i) + dt2 * p(i)
      }
      projectToBase(x)
      projectToTangent(x, p)
    }
  }
}


// solves the equation
//   dx/dt = - dF/dx + sqrt(2 T) eta
//
// using simple Euler integration
//
abstract class OverdampedLangevin(x: Array[R], T: R, dt: R, subIter: Int, rand: Random) {
  val N = x.size
  val f = new Array[R](N)
  
  // fills force array, f = dF/dx, given position x
  def calcForce(x: Array[R], f: Array[R])
  
  // projects vector x onto base manifold
  def projectToBase(x: Array[R]) { }

  def step() {
    val dt2 = dt / subIter
    calcForce(x, f)
    
    for (k <- 0 until subIter) {
      for (i <- 0 until N) {
        x(i) = x(i) - dt2 * f(i) + (math.sqrt(dt2 * 2 * T) * rand.nextGaussian()).toFloat
      }
      projectToBase(x)
    }
    
    if (false) {
      // A way to test derivatives. The expected change in F follows exactly from dF/dS 
      // Note: Noise must be disabled for this to be meaningful.
      val f_dot_f = f.map(x => x*x).sum
          println("expected change in action " + (- f_dot_f * dt2))
    }
  }
}
