package kip.dielectric


import no.uib.cipr.matrix._ 

/**
 * @param order The number of series terms to retain, 'l'.  Order 0 corresponds to pure Coulomb interactions.
 *   Order 1 includes the point-sphere interaction.
 * @param L The center to center distance of the two spheres
 */
// TODO: Place constants into param list
case class Constants(val order: Int, val L: Double) {
  // dielectric constant of background between spheres
  val eps0: Double = 1
  
  // (x_a = eps_a / eps_0) where eps_a is dielectric constant of sphere a
  val (xa, xb) = (10/eps0, 10/eps0)
    
  // (t_a = r_a / L) where r_a is the radius of sphere a
  val (ta, tb) = (5/L, 5/L)
    
  // sphere charges
  val (qa, qb) = (1, -1)
  
  def factorial(i: Int):Double = {
    var ret: Double = 1.0
    for (iter <- 1 to i) ret *= iter
    ret
  }
  def beta(i: Int, j: Int) = factorial(i+j) / (factorial(i)*factorial(j))
  
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

  def forceVector(order: Int, L1: Double, L2: Double) = {
    val dL = (L2 - L1) / 200
    val xs = L1 until L2 by dL
    val fs = xs map (x => force(new Constants(order, x), dL))
    (xs, fs)
  }
  
  def plotAt(order: Int) {
    val (xs, fs) = forceVector(order, 2, 4)
    
    val formatted = kip.lmps.Util.formatDataInColumns(
      ("position", xs.toArray),
      ("force", fs.toArray)
    )
    println(formatted)
  }
  
  def main(args : Array[String]) : Unit = {
    val order = 30
    val Ls = List(11, 12, 13, 15, 17, 20)
//    val Ls = Seq(6.5, 7.2)
    val dL = 1e-4 // seems good if L>=7    
    for (L <- Ls) {
      val c = new Constants(order, L)
      println("L=" + L + " force=" +force(c, 1e-2))
      println("L=" + L + " force=" +force(c, 1e-3))
      println("L=" + L + " force=" +force(c, 1e-4))
    }
  }
}
