package kip.projects.bridging



// A system of equations
//   d_t w + d_x f = 0
// where w and f(w) are vectors representing conserved fields and their fluxes
//
trait ConservationLaws {
  val L: Int
  val ws: Seq[Array[Double]]
  def fs(ws: Seq[Array[Double]]): Seq[Array[Double]]
}


trait CLIntegrator {  
  def halfIter(cl: ConservationLaws, dx: Double, dt: Double)
  
  // Integrate conserved fields cl.ws by timestep dt
  def iter(cl: ConservationLaws, dx: Double, dt: Double) {
    // two half steps on a staggered grid
    halfIter(cl, dx, 0.5*dt)
    halfIter(cl, dx, 0.5*dt)
    
    // rotate back to original index <-> coordinate map
    for (w <- cl.ws) {
      val wn = w(w.size-1)
      for (i <- w.size-1 to 1 by -1) w(i) = w(i-1)
      w(0) = wn
    }
  }
}

object CLIntegratorFirstOrder extends CLIntegrator {  
  def halfIter(cl: ConservationLaws, dx: Double, dt: Double) {
    val fs = cl.fs(cl.ws)
    for ((w, f) <- cl.ws zip fs) {
      val wc = w.clone()
      for (i <- 0 until cl.L) {
        val ip = (i+1) % cl.L
        w(i) = 0.5 * (wc(i) + wc(ip)) - (dt / dx) * (f(ip) - f(i))
      }
    }
  }
}

object CLIntegratorSecondOrder extends CLIntegrator {  
  def minmod(x: Double, y: Double): Double = {
    if (x > 0 == y > 0) 0.5*(x+y)
    else 0
  }
  
  def derivative(w: Array[Double], dx: Double): Array[Double] = {
    val L = w.size
    Array.tabulate(L) { i =>
      val im = (i-1+L)%L
      val ip = (i+1)%L
      minmod(w(ip)-w(i), w(i)-w(im)) / dx 
    }
  }
  
  def halfIter(cl: ConservationLaws, dx: Double, dt: Double) {
    val L = cl.L
    val fs = cl.fs(cl.ws)
    
    // calculate w^(n+1/2)
    val wsp = for ((w, f) <- cl.ws zip fs) yield {
      val df = derivative(f, dx)
      Array.tabulate(L)(i => w(i) - (dt/(2*dx)) * df(i))
    }
    
    // calculate f^(n+1/2)
    val fsp = cl.fs(wsp)
    
    // calculate w^(n+1)
    for (k <- 0 until cl.ws.size) {
      val w  = cl.ws(k)
      val wc = w.clone()
      val dw = derivative(wc, dx)
      val fp = fsp(k)
      
      for (i <- 0 until L) {
        val ip = (i+1)%L
        w(i) = 0.5*(wc(i) + wc(ip)) + (dx/8) * (dw(i) - dw(ip)) - (dt/dx) * (fp(ip) - fp(i))
      }
    }
  }
}


//
// Elastodynamic conservation laws in 1d
//

trait StressFn {
  // Stress as a function of state variables (deformation gradient and energy)
  def stress(defgrad: Double, energy: Double): Double
  
  // Potential energy density at zero temperature
  def zeroTempEnergy(defgrad: Double): Double
  
  // At zero temperature, the stress should be the derivative of the potential energy.
  // Check this by numerical differentiation.
  def checkZeroTempConsistency(defgrad: Double, h: Double) {
    val e = zeroTempEnergy(defgrad)
    val s1 = stress(defgrad, e)
    val s2 = {
      val e1 = zeroTempEnergy(defgrad-h)
      val e2 = zeroTempEnergy(defgrad+h)
      (e2 - e1) / (2*h)
    }
    println("Checking stress at deformation gradient %g. Exact: %g, Numeric: %g".format(defgrad, s1, s2))
  }
}

class ElastodynamicLaws1d(val L: Int, val rho0: Double, defgrad0: Array[Double], sf: StressFn) extends ConservationLaws {
  val defgrad: Array[Double] = defgrad0.clone()
  val momentum: Array[Double] = new Array[Double](L)
  val energy: Array[Double] = new Array[Double](L)
  
  adjustEnergyToZeroTemperature()
  
  // adjust energy so that the configuration is at zero temperature
  def adjustEnergyToZeroTemperature() {
    for (i <- 0 until L) {
      energy(i) = sf.zeroTempEnergy(defgrad(i))
    }
  }
  
  val ws = Seq(defgrad, momentum, energy)
  
  def fs(ws: Seq[Array[Double]]): Seq[Array[Double]] = {
    val w_defgrad = ws(0)
    val w_momentum = ws(1)
    val w_energy = ws(2)
    
    val stress = Array.tabulate(L)(i => sf.stress(defgrad=w_defgrad(i), energy=w_energy(i)))

    val v =     Array.tabulate(L)(i => w_momentum(i) / rho0)
    val sflux = Array.tabulate(L)(i => -v(i))
    val vflux = Array.tabulate(L)(i => -stress(i))
    val eflux = Array.tabulate(L)(i => -stress(i)*v(i))
    
    Seq(sflux, vflux, eflux)
  }
  
}
