package kip.projects.dielectric


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

  def forceVector(consts: Constants, L1: Double, L2: Double) = {
    val dL = (L2 - L1) / 200
    val xs = L1 until L2 by dL
    val fs = xs map (x => force(consts.copy(L=x), dL))
    (xs, fs)
  }
  
  def plotAt(consts: Constants, order: Int) {
    val (xs, fs) = forceVector(consts, 2, 4)
    
    val formatted = kip.util.Util.formatDataInColumns(
      ("position", xs.toArray),
      ("force", fs.toArray)
    )
    println(formatted)
  }
  
  def main(args : Array[String]) : Unit = {
    val c = Constants(
      L     =  0,
      ra    =  5,
      rb    =  5,
      eps0  =  1,
      epsa  = 10,
      epsb  = 10,
      qa    =  1,
      qb    = -1,
      order = 30
    )
    val Ls = List(11, 12, 13, 15, 17, 20)
    val dL = 1e-4 // seems good if L>=7
    
    for (L <- Ls) {
      val cFull = c.copy(L=L)
      val cCoul = c.copy(L=L, epsa=1, epsb=1)
      
      println("L="+L)
      println("Coulomb force  = " +force(cCoul, 1e-3))
      println("Force          = " +force(cFull, 1e-2)+" "+force(cFull, 1e-3)+" "+force(cFull, 1e-4))
      println("Coulomb energy = "+energy(cCoul))
      println("Full energy    = "+energy(cFull))
      println("Delta energy   = "+(energy(cFull) - energy(cCoul)))
      println("Delta fraction = "+(energy(cFull) - energy(cCoul)) / energy(cCoul))
      println()
    }
  }
}
