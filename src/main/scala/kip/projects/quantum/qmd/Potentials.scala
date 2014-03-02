package kip.projects.quantum.qmd

import scala.math.exp
import scala.math.pow
import scala.math.sqrt
import kip.math.Vec3
import Units._

trait Potential {
  def rcut: Double
  
  def phi(r: Double): Double
  def dphi_dr(r: Double): Double
  
  def numOrbitalsPerSite: Int
  
  def fillTBHoppings(del: Vec3,
                     h: Array[Array[Double]],
                     dh_dx: Array[Array[Double]],
                     dh_dy: Array[Array[Double]],
                     dh_dz: Array[Array[Double]])
}

object GoodwinSi extends Potential {
  import math._
  
  // inter-orbital hoppings
  val hssσ = -1.82*eV
  val hspσ = 1.96*eV
  val hpsσ = -hspσ
  val hppσ = 3.06*eV
  val hppπ = -0.87*eV
  
  val en = 2.00           // hopping exponent 
  val em = 4.54           // potential exponent
  
  val Δsp = 8.295*eV      // sp energy separation
  
  val r0 = 2.351*angstrom // equilibrium nearest neighbor distance in Si diamond lattice 
  val rc = 3.67*angstrom  // characteristic truncation length
  val nc = 6.48*angstrom  // sharpness of truncation
  val ϕ1 = 3.4581*eV      // pair potential energy scale  
  
  val rcut = 1.3*rc       // at this distance, hopping decays by factor ~1e-4
  
  
  def phi(r: Double): Double = {
    ϕ1 * pow(r0/r, em) * exp(em * (- pow(r/rc, nc) + pow(r0/rc, nc)))
  }
  
  def dphi_dr(r: Double): Double = {
    - phi(r) * em * (1 + nc*pow(r/rc, nc)) / r
  }
  
  def numOrbitalsPerSite = 4
  
  def fillTBHoppings(del: Vec3,
                     h: Array[Array[Double]],
                     dh_dx: Array[Array[Double]],
                     dh_dy: Array[Array[Double]],
                     dh_dz: Array[Array[Double]]) {
    val x = del.x
    val y = del.y
    val z = del.z
    val r = sqrt(x*x + y*y + z*z)
    
    if (r == 0.0) {
      for (i <- 0 until 4;
           j <- 0 until 4) {
        (if (i != j)
          h(i)(j) = 0.0
        else
          h(i)(j) = if (i == 0) -0.5*Δsp else +0.5*Δsp)
        dh_dx(i)(j) = 0.0
        dh_dy(i)(j) = 0.0
        dh_dz(i)(j) = 0.0
      }
      return
    }
    
    //c<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    //c  compute radial contribution to hopping matrix element
    //c<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    
    val radhop = pow(r0/r, en) * exp(en * (- pow(r/rc, nc) + pow(r0/rc, nc)))
    val dradh = - radhop * en * (1 + nc*pow(r/rc, nc)) / r
    val csssig=hssσ*radhop
    val dsssig=hssσ*dradh
    val cspsig=hspσ*radhop
    val dspsig=hspσ*dradh
    val cppsig=hppσ*radhop
    val dppsig=hppσ*dradh
    val cpppi =hppπ*radhop
    val dpppi =hppπ*dradh
    
    //c<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    //c  compute direction cosines
    //c<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    //c  The minus sign is the phase defined by S-K table in Harrison.
    
    val rinv = 1.0/r
    val xor  = x*rinv
    val l  = -xor
    val l2 = l*l
    val yor  = y*rinv
    val m  = -yor
    val m2 = m*m
    val zor  = z*rinv
    val n  = -zor
    val n2 = n*n
    
    val dldx = -rinv*(m2 + n2)
    val dldy =  rinv*l*m
    val dldz =  rinv*l*n
    val dmdx =  rinv*l*m
    val dmdy = -rinv*(l2 + n2)
    val dmdz =  rinv*m*n
    val dndx =  rinv*l*n
    val dndy =  rinv*m*n
    val dndz = -rinv*(l2 + m2)
    
    val dl2dx = 2.0*l*dldx
    val dl2dy = 2.0*l*dldy
    val dl2dz = 2.0*l*dldz
    val dm2dx = 2.0*m*dmdx
    val dm2dy = 2.0*m*dmdy
    val dm2dz = 2.0*m*dmdz
    val dn2dx = 2.0*n*dndx
    val dn2dy = 2.0*n*dndy
    val dn2dz = 2.0*n*dndz
    
    val dsssx = xor*dsssig
    val dsssy = yor*dsssig
    val dsssz = zor*dsssig
    val dspsx = xor*dspsig
    val dspsy = yor*dspsig
    val dspsz = zor*dspsig
    val dppsx = xor*dppsig
    val dppsy = yor*dppsig
    val dppsz = zor*dppsig
    val dpppx = xor*dpppi
    val dpppy = yor*dpppi
    val dpppz = zor*dpppi
    
    //c<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    //c  s-orbitals with s-orbitals
    //c<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    var csig = 1.0
    var dcsigx = 0.0
    var dcsigy = 0.0
    var dcsigz = 0.0

    var hsig = csig*csssig
    var dhsigx = dcsigx*csssig + csig*dsssx
    var dhsigy = dcsigy*csssig + csig*dsssy
    var dhsigz = dcsigz*csssig + csig*dsssz
    
    var hpi = 0.0
    var dhpix = 0.0
    var dhpiy = 0.0
    var dhpiz = 0.0
    
    h(0)(0) = hsig + hpi
    dh_dx(0)(0) = dhsigx + dhpix
    dh_dy(0)(0) = dhsigy + dhpiy
    dh_dz(0)(0) = dhsigz + dhpiz

    //c<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    //c  s-orbitals with p-orbitals
    //c<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    csig = l
    dcsigx = dldx
    dcsigy = dldy
    dcsigz = dldz

    hsig = csig*cspsig
    dhsigx = dcsigx*cspsig + csig*dspsx
    dhsigy = dcsigy*cspsig + csig*dspsy
    dhsigz = dcsigz*cspsig + csig*dspsz

    hpi = 0.0
    dhpix = 0.0
    dhpiy = 0.0
    dhpiz = 0.0

    h(0)(1) = hsig + hpi
    dh_dx(0)(1) = dhsigx + dhpix
    dh_dy(0)(1) = dhsigy + dhpiy
    dh_dz(0)(1) = dhsigz + dhpiz

    h(1)(0) = - h(0)(1)
    dh_dx(1)(0) =  - dh_dx(0)(1)
    dh_dy(1)(0) =  - dh_dy(0)(1)
    dh_dz(1)(0) =  - dh_dz(0)(1)

    csig = m
    dcsigx = dmdx
    dcsigy = dmdy
    dcsigz = dmdz

    hsig = csig*cspsig
    dhsigx = dcsigx*cspsig + csig*dspsx
    dhsigy = dcsigy*cspsig + csig*dspsy
    dhsigz = dcsigz*cspsig + csig*dspsz

    hpi = 0.0
    dhpix = 0.0
    dhpiy = 0.0
    dhpiz = 0.0

    h(0)(2) = hsig + hpi
    dh_dx(0)(2) = dhsigx + dhpix
    dh_dy(0)(2) = dhsigy + dhpiy
    dh_dz(0)(2) = dhsigz + dhpiz

    h(2)(0) = - h(0)(2)
    dh_dx(2)(0) =  - dh_dx(0)(2)
    dh_dy(2)(0) =  - dh_dy(0)(2)
    dh_dz(2)(0) =  - dh_dz(0)(2)

    csig = n
    dcsigx = dndx
    dcsigy = dndy
    dcsigz = dndz

    hsig = csig*cspsig
    dhsigx = dcsigx*cspsig + csig*dspsx
    dhsigy = dcsigy*cspsig + csig*dspsy
    dhsigz = dcsigz*cspsig + csig*dspsz

    hpi = 0.0
    dhpix = 0.0
    dhpiy = 0.0
    dhpiz = 0.0

    h(0)(3) = hsig + hpi
    dh_dx(0)(3) = dhsigx + dhpix
    dh_dy(0)(3) = dhsigy + dhpiy
    dh_dz(0)(3) = dhsigz + dhpiz

    h(3)(0) = - h(0)(3)
    dh_dx(3)(0) = - dh_dx(0)(3)
    dh_dy(3)(0) = - dh_dy(0)(3)
    dh_dz(3)(0) = - dh_dz(0)(3)
    
    //c<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    //c  p-orbitals with p-orbitals
    //c<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    csig = l2
    dcsigx = dl2dx
    dcsigy = dl2dy
    dcsigz = dl2dz
    
    var cpi = (1.0 - l2)
    var dcpix =  - dl2dx
    var dcpiy =  - dl2dy
    var dcpiz =  - dl2dz
    
    hsig=csig*cppsig
    dhsigx=dcsigx*cppsig + csig*dppsx
    dhsigy=dcsigy*cppsig + csig*dppsy
    dhsigz=dcsigz*cppsig + csig*dppsz

    hpi=cpi*cpppi
    dhpix=dcpix*cpppi + cpi*dpppx
    dhpiy=dcpiy*cpppi + cpi*dpppy
    dhpiz=dcpiz*cpppi + cpi*dpppz

    h(1)(1) = hsig + hpi
    dh_dx(1)(1) = dhsigx + dhpix
    dh_dy(1)(1) = dhsigy + dhpiy
    dh_dz(1)(1) = dhsigz + dhpiz

    csig = l*m
    dcsigx = dldx*m + l*dmdx
    dcsigy = dldy*m + l*dmdy
    dcsigz = dldz*m + l*dmdz

    cpi = -l*m
    dcpix =  - dcsigx
    dcpiy =  - dcsigy
    dcpiz =  - dcsigz

    hsig=csig*cppsig
    dhsigx=dcsigx*cppsig + csig*dppsx
    dhsigy=dcsigy*cppsig + csig*dppsy
    dhsigz=dcsigz*cppsig + csig*dppsz

    hpi=cpi*cpppi
    dhpix=dcpix*cpppi + cpi*dpppx
    dhpiy=dcpiy*cpppi + cpi*dpppy
    dhpiz=dcpiz*cpppi + cpi*dpppz

    h(1)(2) = hsig + hpi
    dh_dx(1)(2) = dhsigx + dhpix
    dh_dy(1)(2) = dhsigy + dhpiy
    dh_dz(1)(2) = dhsigz + dhpiz

    h(2)(1) = hsig + hpi
    dh_dx(2)(1) = dh_dx(1)(2)
    dh_dy(2)(1) = dh_dy(1)(2)
    dh_dz(2)(1) = dh_dz(1)(2)
    
    csig = l*n
    dcsigx = dldx*n + l*dndx
    dcsigy = dldy*n + l*dndy
    dcsigz = dldz*n + l*dndz

    cpi = -l*n
    dcpix = -dcsigx
    dcpiy = -dcsigy
    dcpiz = -dcsigz

    hsig=csig*cppsig
    dhsigx=dcsigx*cppsig + csig*dppsx
    dhsigy=dcsigy*cppsig + csig*dppsy
    dhsigz=dcsigz*cppsig + csig*dppsz

    hpi=cpi*cpppi
    dhpix=dcpix*cpppi + cpi*dpppx
    dhpiy=dcpiy*cpppi + cpi*dpppy
    dhpiz=dcpiz*cpppi + cpi*dpppz

    h(1)(3) = hsig + hpi
    dh_dx(1)(3) = dhsigx + dhpix
    dh_dy(1)(3) = dhsigy + dhpiy
    dh_dz(1)(3) = dhsigz + dhpiz

    h(3)(1) = hsig + hpi
    dh_dx(3)(1) = dh_dx(1)(3)
    dh_dy(3)(1) = dh_dy(1)(3)
    dh_dz(3)(1) = dh_dz(1)(3)

    csig = m2
    dcsigx = dm2dx
    dcsigy = dm2dy
    dcsigz = dm2dz

    cpi = (1.0 - m2)
    dcpix = - dcsigx
    dcpiy = - dcsigy
    dcpiz = - dcsigz

    hsig=csig*cppsig
    dhsigx=dcsigx*cppsig + csig*dppsx
    dhsigy=dcsigy*cppsig + csig*dppsy
    dhsigz=dcsigz*cppsig + csig*dppsz

    hpi=cpi*cpppi
    dhpix=dcpix*cpppi + cpi*dpppx
    dhpiy=dcpiy*cpppi + cpi*dpppy
    dhpiz=dcpiz*cpppi + cpi*dpppz

    h(2)(2) = hsig + hpi
    dh_dx(2)(2) = dhsigx + dhpix
    dh_dy(2)(2) = dhsigy + dhpiy
    dh_dz(2)(2) = dhsigz + dhpiz

    csig = m*n
    dcsigx = dmdx*n + m*dndx
    dcsigy = dmdy*n + m*dndy
    dcsigz = dmdz*n + m*dndz

    cpi = -m*n
    dcpix = - dcsigx
    dcpiy = - dcsigy
    dcpiz = - dcsigz

    hsig=csig*cppsig
    dhsigx=dcsigx*cppsig + csig*dppsx
    dhsigy=dcsigy*cppsig + csig*dppsy
    dhsigz=dcsigz*cppsig + csig*dppsz

    hpi=cpi*cpppi
    dhpix=dcpix*cpppi + cpi*dpppx
    dhpiy=dcpiy*cpppi + cpi*dpppy
    dhpiz=dcpiz*cpppi + cpi*dpppz

    h(2)(3) = hsig + hpi
    dh_dx(2)(3) = dhsigx + dhpix
    dh_dy(2)(3) = dhsigy + dhpiy
    dh_dz(2)(3) = dhsigz + dhpiz

    h(3)(2) = hsig + hpi
    dh_dx(3)(2) = dh_dx(2)(3)
    dh_dy(3)(2) = dh_dy(2)(3)
    dh_dz(3)(2) = dh_dz(2)(3)

    csig = n2
    dcsigx = dn2dx
    dcsigy = dn2dy
    dcsigz = dn2dz

    cpi = (1.0 - n2)
    dcpix = - dcsigx
    dcpiy = - dcsigy
    dcpiz = - dcsigz

    hsig=csig*cppsig
    dhsigx=dcsigx*cppsig + csig*dppsx
    dhsigy=dcsigy*cppsig + csig*dppsy
    dhsigz=dcsigz*cppsig + csig*dppsz

    hpi=cpi*cpppi
    dhpix=dcpix*cpppi + cpi*dpppx
    dhpiy=dcpiy*cpppi + cpi*dpppy
    dhpiz=dcpiz*cpppi + cpi*dpppz

    h(3)(3) = hsig + hpi
    dh_dx(3)(3) = dhsigx + dhpix
    dh_dy(3)(3) = dhsigy + dhpiy
    dh_dz(3)(3) = dhsigz + dhpiz
  }
  
}
