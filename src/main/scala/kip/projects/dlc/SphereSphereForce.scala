package kip.projects.dlc


import no.uib.cipr.matrix._ 


case class Constants(
  order: Int   , // Terms in series. Order 0: pure Coulomb. Order 1: point-sphere correction.
  L:     Double, // Center to center distance between spheres
  ra:    Double, // Radius of sphere a
  rb:    Double, // Radius of sphere b
  eps0:  Double, // Dielectric constant of solution
  epsa:  Double, // Dielectric constant of sphere a
  epsb:  Double, // Dielectric constant of sphere b
  qa:    Double, // Charge of sphere a
  qb:    Double  // Charge of sphere b
  ) {

  val (xa, xb) = (epsa/eps0, epsb/eps0)
  val (ta, tb) = (ra/L, rb/L)
  
  def beta(i: Int, j: Int) = {
    def factorial(i: Int):Double = {
      var ret: Double = 1.0
      for (iter <- 1 to i) ret *= iter
      ret
    }
    factorial(i+j) / (factorial(i)*factorial(j))
  }
  
  def alphaA(l: Double) = (1-xa)*l*math.pow(ta, 2*l+1) / ((xa+1)*l+1)
  def alphaB(l: Double) = (1-xb)*l*math.pow(tb, 2*l+1) / ((xb+1)*l+1)
}


object SphereSphereForce {
  def buildMatrix(consts: Constants) = {
    import consts._
    
    val m = new DenseMatrix(2*order,2*order)
    for (i <- 0 until order; j <- 0 until order) {
      m.set(i,j,              if (i == j) 1 else 0)            // upper left identity matrix
      m.set(order+i, order+j, if (i == j) 1 else 0)            // lower right identity matrix
      m.set(i, order+j,       - alphaA(i+1) * beta(i+1,j+1))   // upper right dense matrix
      m.set(order+i, j,       - alphaB(i+1) * beta(i+1,j+1))   // lower left dense matrix
    }
    m
  }
  
  def buildRHSVector(consts: Constants) = {
    import consts._
    
    val v = new DenseVector(2*order)
    for (i <- 0 until order) {
      v.set(i,       qb * alphaA(i+1) / eps0)
      v.set(i+order, qa * alphaB(i+1) / eps0)
    }
    v
  }
  
  def energy(consts: Constants) = {
    import consts._
    
    val m = buildMatrix(consts)
    val v = buildRHSVector(consts)
    val qvec = new DenseVector(v.size)
    m.solve(v, qvec)
    
    var energy = qa * qb / (eps0*L)
    for (i <- 0 until order) {
      energy += (1/(2*L)) * (qb*qvec.get(i) + qa*qvec.get(i+order))
    }
    
    val correction = 1 / (4*math.Pi) // switch units
    energy * correction
  }
  
  def force(consts: Constants, dL: Double) = {
    val e1 = energy(consts.copy(L=consts.L-dL/2))
    val e2 = energy(consts.copy(L=consts.L+dL/2))
    - (e2-e1) / dL
  }


  def makePlot() {
    // Two spheres
    // val c = Constants(
    //   L     =  0,
    //   ra    =  1,
    //   rb    =  1,
    //   eps0  =  1,
    //   epsa  = 10,
    //   epsb  = 10,
    //   qa    =  1,
    //   qb    = -1,
    //   order = 30
    // )
    // val L1 = 2.001
    // val L2 = 2.200

    // Messina's system
    val c = Constants(
      L     =  0,
      ra    =  7.5,
      rb    =  0,
      eps0  =  80,
      epsa  =  2,
      epsb  =  0,
      qa    = -60,
      qb    =  3,
      order = 40
    )
    val L1 = 8.00
    val L2 = 15
    
    val dL = 0.05
    val orders = Seq(45, 60, 100)
    
    val xs = (L1 to L2 by dL).toArray
    val us = orders map (o => xs map (x => energy(c.copy(L=x, order=o))))
    val fs = orders map (o => xs map (x => force (c.copy(L=x, order=o), dL)))
    val formatted = kip.util.Util.formatDataInColumns(
      ("position", xs),
      ("energy1", us(0)),
      ("energy2", us(1)),
      ("energy3", us(2)),
      ("force1", fs(0)),
      ("force2", fs(1)),
      ("force3", fs(2)),
      ("energyCoul", xs.map (x => energy(c.copy(L=x, order=0))))
    )
    
    println(formatted)
    // kip.util.Util.writeStringToFile(formatted, "yahoo.txt")
  }
  
  def printFew() {
    // Daniel's system
    val c = Constants(
      L     =  0,
      ra    =  5,
      rb    =  0,
      eps0  =  1,
      epsa  = 10,
      epsb  =  1,
      qa    =  0,
      qb    =  1,
      order = 40
    )
    val Ls = List(6, 7, 8, 9, 10, 12, 15)
    
    for (L <- Ls) {
      val cFull = c.copy(L=L)
      val cCoul = c.copy(L=L, epsa=1, epsb=1)

      println("L="+L)
      println("Coulomb force  = %e" format force(cCoul, 1e-3))
      println("Force          = %e" format force(cFull, 1e-3))
      println("Coulomb energy = %e" format energy(cCoul))
      println("Full energy    = %e" format energy(cFull))
      println("Delta energy   = %e" format (energy(cFull) - energy(cCoul)))
      println("Delta fraction = %e" format (energy(cFull) - energy(cCoul)) / energy(cCoul))
      println()
    }
  }
  
  def smallLargePRL() {
    val c0 = Constants(
        L     =  4,
        ra    =  3,
        rb    =  0.5,
        eps0  =  1,
        epsa  =  1,
        epsb  =  1,
        qa    =  1,
        qb    = -1,
        order = 40
        )
    val c1 = c0.copy(epsa=100, epsb=1)
    val c2 = c0.copy(epsa=100, epsb=100)
    val c = Seq(c0, c1, c2)
    val e = c.map(4*math.Pi*SphereSphereForce.energy(_))
    val f = c.map(4*math.Pi*SphereSphereForce.force(_, 1e-3))
    
    println(s"q_a=${c0.qa}, q_b=${c0.qb}")
    println()
    println("kappa_a=1, kappa_b=1")
    println("  energy0="+ e(0))       // -1/4
    println("  force0 ="+ f(0))       // -1/16
    
    println("kappa_a=100, kappa_b=1")
    println("  energy1="+ e(1))       // -0.367403
    println("  force1 ="+ f(1))       // -0.255647

    println("kappa_a=100, kappa_b=100")
    println("  energy2="+ e(2))       // -0.371822
    println("  force2 ="+ f(2))       // -0.273539
    
    println("relative effect of kappa_b")
    println("  (E2-E1)/E2 = " + (e(2) - e(1)) / e(2))
    println("  (f2-f1)/f2 = " + (f(2) - f(1)) / f(2))
  }
  
  def main(args : Array[String]) : Unit = {
    // printFew()
    //makePlot()
    smallLargePRL()
  }
}
