package kip.projects.quantum

import smatrix._

object ScalarTyp {
  type S = Scalar.ComplexDbl
  type R = S#Raw
  val ctor = Constructors.complexDbl
}
import ScalarTyp._
import ctor._


object Quantum extends App {
  // testIntegratedDensity()
  testDerivative()
  
  // Calculates effective action at given filling fraction for various configurations
  def testEigenvalues() {
    val q = new Quantum(w=20, h=20, t=1, J_eff=0.3, e_min= -10, e_max=10)
    val n = q.matrix.numRows
    println("Matrix dim = "+n)
    
    q.setFieldAllOut(q.field)
    q.fillMatrix(q.matrix)
    val eig1 = KPM.eigenvaluesExact(q.matrix)
    val i_cut = n * 3 / 4
    val e_cut = eig1(i_cut-1) // (eig1(i_cut-1) + eig1(i_cut)) / 2.0
    val filledEig1 = eig1.takeWhile(_ < e_cut)
    println("Gap between: [%g, %g]".format(eig1(i_cut-1), eig1(i_cut)))
    println("Choosing potential mu = %g".format(e_cut))
    println("All-out fraction = %g".format(filledEig1.size / n.toDouble))
    println("All-out weight = %g".format(filledEig1.map(_.re - e_cut).sum))
    
    q.setFieldFerro(q.field)
    q.fillMatrix(q.matrix)
    val eig2 = KPM.eigenvaluesExact(q.matrix)
    val filledEig2 = eig2.takeWhile(_ < e_cut)
    println("Ferro fraction = %g".format(filledEig2.size / n.toDouble))
    println("Ferro weight = %g".format(filledEig2.map(_.re - e_cut).sum))
  }
  
  def testSpeed {
    val q = new Quantum(w=10, h=10, t=1, J_eff=2, e_min= -10, e_max= 10)  // hopping only: e_min= -6-0.5, e_max= 3+0.5
    val H = q.matrix
    val kpm = new KPM(H, order=100, nrand=1)
    
    val r = kpm.randomVector()
    val c = kpm.expansionCoefficients(de=1e-4, e => e)
 
    val dH = H.duplicate
    for (i <- 0 until 10) {
      kip.util.Util.time("Forward")(kpm.momentsStochastic(r))
      kip.util.Util.time("Backward")(kpm.functionAndGradient(r, c, dH))
    }
  }

  // Plots the integrated density of states
  def testIntegratedDensity() {
    val q = new Quantum(w=10, h=10, t=1, J_eff=2, e_min= -10, e_max= 10)  // hopping only: e_min= -6-0.5, e_max= 3+0.5
    val H = q.matrix
    require((H - H.dag).norm2.abs < 1e-10, "Found non-hermitian hamiltonian!")
    println("N = "+H.numRows)
    val kpm = new KPM(H, order=100, nrand=1)
    val range = kpm.range

      val plot = KPM.mkPlot("Integrated density of states")
      KPM.plotLines(plot, (kpm.range, KPM.integrateDeltas(range, KPM.eigenvaluesExact(H), moment=0)), "Exact", java.awt.Color.RED)
      KPM.plotLines(plot, (kpm.range, KPM.integrate(range, kpm.eigenvaluesApprox(kpm.jacksonKernel), moment=0)), "Approx", java.awt.Color.BLACK)
  }
  
  def testDerivative() {
    val q = new Quantum(w=10, h=10, t=1, J_eff=2, e_min= -10, e_max= 10)  // hopping only: e_min= -6-0.5, e_max= 3+0.5
    val H = q.matrix
    val dH = q.delMatrix
    val kpm = new KPM(H, order=100, nrand=1)
    val r = kpm.randomVector()
    val c = kpm.expansionCoefficients(de=1e-4, e => e*e)
    val f0 = kpm.functionAndGradient(r, c, dH)
    println("H = "+H)
    println("dH = "+dH)

    val k = 13
    val del = 1e-7
    q.fieldDerivative(dH, q.delField)
    val deriv = q.delField(k)
    
    q.field(k) += del
    q.fillMatrix(H)
    println("new H = "+H)
    
    val f1 = kpm.functionAndGradient(r, c, dH)
    println("deriv = "+deriv)
    println("raw fn: f0 = %g f1 = %g".format(f0, f1))
    println("approx deriv: (f1 - f0)/del = "+ (f1 - f0)/del)
    println("error1: (f1 - f0)/del - dH = "+((f1 - f0)/del - deriv))
  }
}


// Notation:
//  d    = vector component (3 dimensional)
//  sp   = Dirac spin index
//  x, y = coordinates on triangular lattice
//  n    = nearest neighbor index on lattice 
//
class Quantum(val w: Int, val h: Int, val t: R, val J_eff: R, val e_min: R, val e_max: R) {
  require(h % 2 == 0, "Need even number of rows, or hamiltonian is non-hermitian")
  val vectorDim = 3
  
  val e_avg   = (e_max + e_min)/2
  val e_scale = (e_max - e_min)/2
  
  def matrixIndex(sp: Int, x: Int, y: Int): Int = {
    sp + x*(2) + y*(2*w)
  }
  
  def fieldIndex(d: Int, x: Int, y: Int): Int = {
    d + x*(3) + y*(3*w)
  }

  def pauliIndex(sp1: Int, sp2: Int, d: Int): Int = {
    sp1 + sp2*(2) + d*(2*2)
  }
  def pauli = Array[S#A] (
    0, 1,
    1, 0,
    
    0, I, // visually transposed, due to row major ordering
   -I, 0,
    
    1, 0,
    0, -1
  )
  
  val field: Array[R] = {
    val ret = new Array[R](vectorDim*w*h)
    setFieldFerro(ret)
    ret
  }
  val delField: Array[R] = Array.fill(vectorDim*w*h)(0)
  
  val matrix = {
    val ret = sparse(2*h*w, 2*h*w)
    fillMatrix(ret)
    ret.toPacked
  }
  val delMatrix: PackedSparse[S] = {
    matrix.duplicate.clear
  }
  
  def setFieldFerro(field: Array[R]) { 
    field.transform(_ => 1.0)
    normalizeField(field)
  }
  
  def setFieldRandom(field: Array[R], rand: util.Random) {
    field.transform(_ => rand.nextDouble() - 0.5)
    normalizeField(field)
  }
  
  def setFieldAllOut(field: Array[R]) {
    for (x <- 0 until w;
         y <- 0 until h) {
      val s = (x%2, y%2) match {
        case (0, 0) => Seq(+1, +1, +1)
        case (0, 1) => Seq(+1, -1, -1)
        case (1, 0) => Seq(-1, +1, -1)
        case (1, 1) => Seq(-1, -1, +1)
      }
      for (d <- 0 until vectorDim) { 
        field(fieldIndex(d, x, y)) = s(d)
      }
    }
    normalizeField(field)
  }
  
  def normalizeField(field: Array[R], validate: Boolean = false) {
    for (y <- 0 until h;
         x <- 0 until w) {
      var acc = 0d
      for (d <- 0 until 3) {
        acc += field(fieldIndex(d, x, y)).abs2
      }
      acc = math.sqrt(acc)
      if (validate)
        require(acc > 0.5 && acc < 1.5, "Vector deviates too far from normalization")
      for (d <- 0 until 3) {
        field(fieldIndex(d, x, y)) /= acc
      }
    }
  }
  
  //   
  //   o - o - o
  //   | \ | \ |    y
  //   o - o - o    ^
  //   | / | / |    |
  //   o - o - o    ----> x
  //
  def neighbors(x: Int, y: Int): (Array[Int], Array[Int]) = {
    val xdel = (y % 2) match {

      //      2   1
      //      | /
      //  3 - o - 0
      //      | \
      //      4   5
      case 0 => Seq(1, 1, 0, -1, 0, 1)

      //  2   1
      //    \ |
      //  3 - o - 0
      //    / |
      //  4   5
      case 1 => Seq(1, 0, -1, -1, -1, 0)
    }
    val ydel = Seq(0, 1, 1, 0, -1, -1)
    
    val xs = Array.tabulate(6) { d => (x+xdel(d)+w)%w }
    val ys = Array.tabulate(6) { d => (y+ydel(d)+h)%h } 
    (xs, ys)
  }

    
  //   
  //   o - o - o
  //   | / | / |    y
  //   o - o - o    ^
  //   | / | / |    |
  //   o - o - o    ----> x
  //
  def neighbors2(x: Int, y: Int): (Array[Int], Array[Int]) = {
    //      2   1
    //      | /
    //  3 - o - 0
    //    / |
    //  4   5
    val xdel = Seq(1, 1, 0, -1, -1, 0)
    val ydel = Seq(0, 1, 1, 0, -1, -1)
    val xs = Array.tabulate(6) { d => (x+xdel(d)+w)%w }
    val ys = Array.tabulate(6) { d => (y+ydel(d)+h)%h } 
    (xs, ys)
  }
  
  def fillMatrix[M[s <: Scalar] <: Sparse[s, M]](m: M[S]) {
    m.clear()
    
    // loop over all lattice sites
    for (y <- 0 until h;
         x <- 0 until w) {
      
      // neighbor hopping
      val (xn, yn) = neighbors2(x, y)
      for (n <- 0 until 6;
           sp <- 0 until 2) {
        val i = matrixIndex(sp, x, y)
        val j = matrixIndex(sp, xn(n), yn(n))
        m(i, j) = -t
      }
      
      // hund coupling 
      for (sp1 <- 0 until 2;
           sp2 <- 0 until 2) {
        
        var coupling = 0: S#A
        for (d <- 0 until 3) {
          coupling += pauli(pauliIndex(sp1, sp2, d)) * field(fieldIndex(d, x, y))
        }
        val i = matrixIndex(sp1, x, y)
        val j = matrixIndex(sp2, x, y)
        m(i, j) = -J_eff * coupling
      }
    }
    
    // scale matrix appropriately so that eigenvalues lie beween -1 and +1
    for (i <- 0 until m.numRows) { m(i,i) -= e_avg }
    m /= e_scale
  }

  def scaleEnergy(x: R): R = {
    (x - e_avg) / e_scale
  }
  
  // Derivative of field S corresponding to derivative of matrix H (hund coupling)
  def fieldDerivative(dH: PackedSparse[S], dS: Array[R]) {
    // loop over all lattice sites and vector indices
    for (y <- 0 until h;
         x <- 0 until w;
         d <- 0 until 3) {
      
      var dCoupling: S#A = 0
      for (sp1 <- 0 until 2;
           sp2 <- 0 until 2) {
        val i = matrixIndex(sp1, x, y)
        val j = matrixIndex(sp2, x, y)
        dCoupling += dH(i, j) * pauli(pauliIndex(sp1, sp2, d))
      }
      require(math.abs(dCoupling.im) < 1e-10, "Imaginary part of field derivative non-zero")
      dS(fieldIndex(d, x, y)) = -J_eff * dCoupling.re
    }
    
    // consistent matrix scaling
    dS.transform(_ / e_scale)
  }
}
