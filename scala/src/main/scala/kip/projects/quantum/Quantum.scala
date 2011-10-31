package kip.projects.quantum

import smatrix._
import ctor._


object Quantum extends App {
  testIntegratedDensity()
//  testDerivative2()
//  testEigenvalues()
  
  // Calculates effective action at given filling fraction for various configurations
  def testEigenvalues() {
    val q = new Quantum(w=30, h=30, t=1, J_H=0.2, e_min= -10, e_max=10)
    val n = q.matrix.numRows
    println("Matrix dim = "+n)
    
    q.setFieldAllOut(q.field)
    q.fillMatrix(q.matrix)
    var eig = KPM.eigenvaluesExact(q.matrix)
    val i_cut = n * 3 / 4
    val e_cut1 = eig(i_cut-1)
    val e_cut2 = (eig(i_cut-1) + eig(i_cut)) / 2.0
    println("Gap between: [%g, %g]".format(eig(i_cut-1), eig(i_cut)))
    println("Low = %g, mid = %g".format(e_cut1, e_cut2))
    println()
    
    def weight(eig: Array[Double], cut: Double) = eig.takeWhile(_ <= cut).map(_.re - cut).sum
    def filling(eig: Array[Double], cut: Double) = eig.takeWhile(_ <= cut).size / n.toDouble
    
    println("All-out   low = %g  mid %g  (frac: %g, %g)".format(weight(eig, e_cut1), weight(eig, e_cut2),filling(eig, e_cut1), filling(eig, e_cut2)))
    
    q.setFieldThreeOut(q.field)
    q.fillMatrix(q.matrix)
    eig = KPM.eigenvaluesExact(q.matrix)
    println("Three-out low = %g  mid %g  (frac: %g, %g)".format(weight(eig, e_cut1), weight(eig, e_cut2),filling(eig, e_cut1), filling(eig, e_cut2)))

    q.setFieldFerro(q.field)
    q.fillMatrix(q.matrix)
    eig = KPM.eigenvaluesExact(q.matrix)
    println("Ferro     low = %g  mid %g  (frac: %g, %g)".format(weight(eig, e_cut1), weight(eig, e_cut2),filling(eig, e_cut1), filling(eig, e_cut2)))
  }
  
  def testSpeed {
    val q = new Quantum(w=10, h=10, t=1, J_H=2, e_min= -10, e_max= 10)
    val H = q.matrix
    val kpm = new KPM(H, nrand=1)
    val order = 100
    val r = kpm.randomVector()
    val c = kpm.expansionCoefficients(order, de=1e-4, e => e)
 
    val dH = H.duplicate
    for (i <- 0 until 10) {
      kip.util.Util.time("Forward")(kpm.momentsStochastic(order, r))
      kip.util.Util.time("Backward")(kpm.functionAndGradient(r, c, dH))
    }
  }

  // Plots the integrated density of states
  def testIntegratedDensity() {
    val q = new Quantum(w=8, h=8, t=1, J_H=0.0, e_min= -10, e_max= 10)
    val H = q.matrix
    require((H - H.dag).norm2.abs < 1e-10, "Found non-hermitian hamiltonian!")
    println("N = "+H.numRows)
    val order = 100
    val kpm = new KPM(H, nrand=1)
    val range = KPM.range(npts=5*order)
    
    val plot = KPM.mkPlot("Integrated density of states")
    KPM.plotLines(plot, (range, KPM.integrateDeltas(range, KPM.eigenvaluesExact(H), moment=0)), "Exact", java.awt.Color.RED)
    KPM.plotLines(plot, (range, KPM.integrate(range, KPM.eigenvaluesApprox(order, kpm), moment=0)), "Approx", java.awt.Color.BLACK)
  }
  
  def testDerivative() {
    val q = new Quantum(w=10, h=10, t=1, J_H=2, e_min= -10, e_max= 10)
    val H = q.matrix
    val dH = q.delMatrix
    val order = 100
    val kpm = new KPM(H, nrand=1)
    val r = kpm.randomVector()
    val c = kpm.expansionCoefficients(order, de=1e-4, e => e*e)
    val f0 = kpm.functionAndGradient(r, c, dH)
    println("H = "+H)
    println("dH = "+dH)

    val k = 13
    val del = 1e-7
    q.fieldDerivative(dH, q.delField)
    val deriv = q.delField(k)
    
    q.field(k) += del
    q.normalizeField(q.field)
    q.fillMatrix(H)
    println("new H = "+H)
    
    val f1 = kpm.functionAndGradient(r, c, dH)
    println("deriv = "+deriv)
    println("raw fn: f0 = %g f1 = %g".format(f0, f1))
    println("approx deriv: (f1 - f0)/del = "+ (f1 - f0)/del)
    println("error1: (f1 - f0)/del - dH = "+((f1 - f0)/del - deriv))
  }
  
  
  def testDerivative2() {
    val q = new Quantum(w=10, h=10, t=1, J_H=0.1, e_min= -10, e_max= 10)
    val H = q.matrix
    val dH1 = q.delMatrix
    val dH2 = dH1.duplicate
    val order = 100
    val kpm = new KPM(H, nrand=200)
    
    val c = kpm.expansionCoefficients(order, de=1e-4, e => e*e)
    kpm.gradientExact(c, dH1)
    val r = kpm.randomGaussianVector()
    kpm.functionAndGradient(r, c, dH2)
    
    println("exact  dH = "+dH1)
    println("approx dH = "+dH2)
    
    val resid = dH1 - dH2
    println("scaled error = " + math.sqrt(resid.norm2.re / dH1.norm2.re))
  }
}


// Notation:
//  d    = vector component (3 dimensional)
//  sp   = Dirac spin index
//  x, y = coordinates on triangular lattice
//  n    = nearest neighbor index on lattice 
//
class Quantum(val w: Int, val h: Int, val t: R, val J_H: R, val e_min: R, val e_max: R) {
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
  
  def setFieldThreeOut(field: Array[R]) {
    for (x <- 0 until w;
         y <- 0 until h) {
      val s = (x%2, y%2) match {
        case (0, 0) => Seq(-1, -1, -1)
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
      if (validate && !(acc > 0.95 && acc < 1.05))
        println("Vector magnitude %g deviates too far from normalization".format(acc))
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
        m(i, j) = -J_H * coupling
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
      require(math.abs(dCoupling.im) < 1e-5, "Imaginary part of field derivative non-zero")
      dS(fieldIndex(d, x, y)) = -J_H * dCoupling.re
    }
    
    // remove component of dS parallel to S
    for (y <- 0 until h;
         x <- 0 until w) {
      var s_dot_s = 0d
      var s_dot_ds = 0d
      for (d <- 0 until 3) {
        val i = fieldIndex(d, x, y)
        s_dot_s  += field(i)*field(i)
        s_dot_ds += field(i)*dS(i) 
      }
      val alpha = s_dot_ds / s_dot_s
      for (d <- 0 until 3) {
        val i = fieldIndex(d, x, y)
        dS(i) -= alpha * field(i)
      }
    }
    
    // consistent matrix scaling
    dS.transform(_ / e_scale)
  }
}
