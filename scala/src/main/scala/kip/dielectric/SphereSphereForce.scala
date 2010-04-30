package kip.dielectric


// This class requires Scalala, which is a pain to compile.

/*
import scalala.Scalala._;
import scalala.tensor.{Tensor,Tensor1,Tensor2,Vector,Matrix};
import scalala.tensor.dense.{DenseVector,DenseMatrix};
import scalala.tensor.sparse._;

case class Constants(val order: Int, val L: Double) {
  val eps0: Double = 1
  
  val (xa, xb) = (10/eps0, 10/eps0)
  val (ta, tb) = (5/L, 5/L)
  val (qa, qb) = (1, 1)
  
  def factorial(i: Int):Double = {
    var ret: Double = 1.0
    for (iter <- 1 to i) ret *= iter
    ret
  }
  def beta(i: Int, j: Int) = factorial(i+j) / (factorial(i)*factorial(j))
  
  def alphaA(l: Double) = (1-xa)*l*pow(ta, 2*l+1) / ((xa+1)*l+1)
  def alphaB(l: Double) = (1-xb)*l*pow(tb, 2*l+1) / ((xb+1)*l+1)
}


object SphereSphereForce {
  def buildMatrix(consts: Constants) = {
    import consts._    
    
    val m = DenseMatrix(2*order,2*order)(0.0)
    for (i <- 0 until order; j <- 0 until order) {
      m(i,j) = if (i == j) 1 else 0            // upper left identity matrix
      m(order+i, order+j) = if (i == j) 1 else 0      // lower right identity matrix
      m(i, order+j) = - alphaA(i+1) * beta(i+1,j+1)  // upper right dense matrix
      m(order+i, j) = - alphaB(i+1) * beta(i+1,j+1)  // lower left dense matrix
    }
    m
  }
  
  def buildRHSVector(consts: Constants) = {
    import consts._
    
    val v = DenseVector(2*order)(0.0)
    for (i <- 0 until order) {
      v(i)    = qb * alphaA(i+1) / eps0
      v(i+order)  = qa * alphaB(i+1) / eps0 
    }
    v
  }
  
  def energy(consts: Constants) = {
    import consts._
    
    val m = buildMatrix(consts)
    val v = buildRHSVector(consts)
    val qvec = m \ v
    
    var energy = qa * qb / (eps0*L)
    for (i <- 0 until order) {
      energy += (1/(2*L)) * (qb*qvec(i) + qa*qvec(i+order))
    }
    
    val correction = 1 / (4*java.lang.Math.PI) // switch units
    energy * correction
  }
  
  def force(consts: Constants, dL: Double) = {
    val e1 = energy(consts.copy(L=consts.L-dL/2))
    val e2 = energy(consts.copy(L=consts.L+dL/2))
    - (e2-e1) / dL
  }

  def forceVector(order: Int, L1: Double, L2: Double) = {
    val x = linspace(L1, L2)
    val y = x.copy
    for (i <- 0 until x.size) {
      y(i) = force(new Constants(order, x(i)), 1e-3)
    }
    (x, y)
  }
  
  def main(args : Array[String]) : Unit = {
//    val orders = List(0, 5, 10, 20, 30)
//    for (order <- orders) {
//      val (x, y) = forceVector(order, 11, 20)
//      scikit.util.Commands.replot(new PointSet(x.toArray, y.toArray), order.toString)
//    }
    
    val order = 30
    val Ls = List(11, 12, 13, 15, 17, 20)
    val dL = 1e-4 // seems good if L>=7    
    for (L <- Ls) {
      val c = new Constants(order, L)
      println("L=" + L + " force=" +force(c, 1e-2))
      println("L=" + L + " force=" +force(c, 1e-3))
      println("L=" + L + " force=" +force(c, 1e-4))
    }
  }
  
  def plotAt(order: Int) {
    val (x, y) = forceVector(order, 2, 4)
    plot(x, y)
  }
}

*/
