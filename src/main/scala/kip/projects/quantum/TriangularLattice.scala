package kip.projects.quantum

import smatrix._
import ctor._


object TriangularLattice extends App {
  import kip.util.Util.time

//  time("eigenvalues")(testEigenvalues())
//  time("speed")(testSpeed())
//  time("integrated density")(testIntegratedDensity())
//  time("derivative")(testDerivative())
  testEffectiveTemperature()
  
  
  // Calculates effective action at given filling fraction for various configurations
  def testEigenvalues() {
//    val q = new TriangularLattice(w=8, h=8, t=1, J_H=2.0, B_n= 0, e_min= -10, e_max=10)
    val q = new TriangularLattice(w=20, h=20, t=1, J_H=3, B_n= 1, e_min= -10, e_max=10)
    val n = q.matrix.numRows
    println("Matrix dim = "+n)
    
    q.setFieldAllOut(q.field)
    q.fillMatrix(q.matrix)
    var eig = KPM.eigenvaluesExact(q.matrix)
    val i_cut = n * 1 / 4
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
    val q = new TriangularLattice(w=10, h=10, t=1, J_H=2, B_n=0, e_min= -10, e_max= 10)
    val H = q.matrix
    val kpm = new KPM(H, nrand=1)
    val order = 100
    val r = kpm.randomVector()
    val c = KPM.expansionCoefficients2(order, quadPts=10*order, e => e)
 
    val dH = H.duplicate
    for (i <- 0 until 10) {
      kip.util.Util.time("Forward")(kpm.momentsStochastic(order, r))
      kip.util.Util.time("Backward")(kpm.functionAndGradient(r, c, dH))
    }
  }

  // Plots the integrated density of states
  def testIntegratedDensity() {
    val q = new TriangularLattice(w=20, h=20, t=1, J_H=3, B_n= 1, e_min= -10, e_max= 10)
    q.setFieldAllOut(q.field)
    q.fillMatrix(q.matrix)
    
    val H = q.matrix
    require((H - H.dag).norm2.abs < 1e-10, "Found non-hermitian hamiltonian!")
    println("N = "+H.numRows)
    val order = 100
    val kpm = new KPM(H, nrand=1, seed=0)
    val range = KPM.range(npts=5*order)
    
    val plot = KPM.mkPlot("Integrated density of states")
    KPM.plotLines(plot, (range, KPM.integrateDeltas(range, KPM.eigenvaluesExact(H), moment=0)), "Exact", java.awt.Color.RED)
    KPM.plotLines(plot, (range, KPM.integrate(range, KPM.eigenvaluesApprox(order, range, kpm), moment=0)), "Approx", java.awt.Color.BLACK)
  }
  
  def testDerivative() {
    val q = new TriangularLattice(w=10, h=10, t=1, J_H=2, B_n=0 ,e_min= -10, e_max= 10)
    val H = q.matrix
    val dH = q.delMatrix
    val order = 100
    val kpm = new KPM(H, nrand=1)
    val r = kpm.randomVector()
    val c = KPM.expansionCoefficients2(order, quadPts=10*order, e => e*e)
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
  
  // Estimate effective temperature as a function of dt_per_rand and J_H
  // (Other parameters shouldn't have much effect)
  // Note: for small J, we observe (T_eff ~ (0.03 J)^2 dt_per_rand) in scaled units.
  // In unscaled units: T => 10 T, dt => 0.1 dt, so that (T_eff ~ (0.3 J)^2 dt)
  // 
  def testEffectiveTemperature() {
    val dt_per_rand = 0.1
    val nrand = 4
    val dt = dt_per_rand * nrand
    val mu = 0 // -2.56
    val q = new TriangularLattice(w=16, h=16, t=1, J_H=2.0, B_n=0, e_min= -10, e_max= 10)
    q.setFieldRandom(q.field, new util.Random(seed=0))
//    q.setFieldFerro(q.field)
    q.fillMatrix(q.matrix)
    
    val H = q.matrix
    val dH1 = q.delMatrix
    val dH2 = dH1.duplicate
    val order = 100
    val kpm = new KPM(H, nrand, seed=2)
    
    val mup = q.scaleEnergy(mu)
    val c = KPM.expansionCoefficients2(order, quadPts=10*order, (e => if (e < mup) (e - mup) else 0))
    kpm.gradientExactDense(c, dH1)
    val r = kpm.randomGaussianVector()
    kpm.functionAndGradient(r, c, dH2)
    
    val dS1 = q.delField
    val dS2 = dS1.clone()
    
    println("exact  dH = "+dH1)
    println("approx dH = "+dH2)

    q.fieldDerivative(dH1, dS1)
    q.fieldDerivative(dH2, dS2)
    
    // get deviations from exact gradient components
    val delta = (dS1, dS2).zipped.map((s1, s2) => s1-s2)
    // sigma = standard deviation
    val sigma2 = delta.map(x => x*x).sum / delta.size
    val sigma = math.sqrt(sigma2)
    
    println("Standard deviation of energy gradient = " + sigma)
    println("Effective (scaled) temperature = " + (dt * sigma2 / 2))
    
    println("""Expected output
Standard deviation of energy gradient = 0.03606297374255914
Effective (scaled) temperature = 2.6010761503130196E-4
""")

    if (false) {
      // plot histogram of error, and compare to gaussian
      import math._
      val hist = new scikit.dataset.Histogram(0.04)
      hist.setNormalizing(true)
      delta.foreach(e => hist.accum(e))
      val gauss = new scikit.dataset.Function { def eval(x: Double) = (1/(sigma*sqrt(2*Pi)))*exp(-x*x/(2*sigma2)) }
      scikit.util.Commands.plot(hist)
      scikit.util.Commands.replot(gauss)
    }
  }
}



// Notation:
//  d    = vector component (3 dimensional)
//  sp   = Dirac spin index
//  x, y = coordinates on triangular lattice
//  nn   = [0,1,...5] nearest neighbor index (oriented clockwise, starting at 3 o'clock)
//
class TriangularLattice(val w: Int, val h: Int, val t: R, val J_H: R, val B_n: Int, val e_min: R, val e_max: R) extends KondoHamiltonian {
  require(h % 2 == 0, "Need even number of rows, or hamiltonian is non-hermitian")

  val numLatticeSites = w*h
  
  def latticeIndex(x: Int, y: Int): Int = {
    x + y*(w)
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
      for (d <- 0 until 3) { 
        field(fieldIndex(d, latticeIndex(x, y))) = s(d)
      }
    }
    normalizeField(field)
  }
  
  def setField1q(field: Array[R]) {
    for (x <- 0 until w;
         y <- 0 until h) {
      val s = (x%2, y%2) match {
        case (0, 0) => Seq(-1, 0, 0)
        case (0, 1) => Seq(-1, 0, 0)
        case (1, 0) => Seq(+1, 0, 0)
        case (1, 1) => Seq(+1, 0, 0)
      }
      for (d <- 0 until 3) { 
        field(fieldIndex(d, latticeIndex(x, y))) = s(d)
      }
    }
    normalizeField(field)
  }
  
  def setField2q(field: Array[R]) {
    for (x <- 0 until w;
         y <- 0 until h) {
      val s = (x%2, y%2) match {
        case (0, 0) => Seq(-1, -1, 0)
        case (0, 1) => Seq(-1, +1, 0)
        case (1, 0) => Seq(+1, -1, 0)
        case (1, 1) => Seq(+1, +1, 0)
      }
      for (d <- 0 until 3) { 
        field(fieldIndex(d, latticeIndex(x, y))) = s(d)
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
      for (d <- 0 until 3) { 
        field(fieldIndex(d, latticeIndex(x, y))) = s(d)
      }
    }
    normalizeField(field)
  }

    //
    //  o - o - o
    //   \ / \ / \    y
    //    o - o - o    ^
    //     \ / \ / \    \
    //      o - o - o    ----> x
    //
   def position(x: Int, y: Int): (Double, Double) = {
    val a = 1                                   // horizontal distance between columns
    val b = 0.5*math.sqrt(3.0)*a                // vertical distance between rows
	(a*x - 0.5*a*y, b*y)
  }
  
  // returns (x, y) indices for the `nn`th neighbor site
  def neighbors(x: Int, y: Int, nn: Int): (Int, Int) = {
    //      2   1
    //      | /
    //  3 - o - 0
    //    / |
    //  4   5
    val xdel = Seq(1, 1, 0, -1, -1, 0)
    val ydel = Seq(0, 1, 1, 0, -1, -1)
    ((x+xdel(nn)+w)%w, (y+ydel(nn)+h)%h)
  }
  
  // returns (dx, dy) positions
  def displacement(x: Int, y: Int, nn: Int): (Double, Double) = {
    val a = 0.5*math.sqrt(3.0)
    val xdisp = Seq(1, 0.5, -0.5, -1, -0.5, 0.5)
    val ydisp = Seq(0, a, a, 0, -a, -a)
    (xdisp(nn), ydisp(nn))
  }
  
  // in *unscaled* (physical) energy units 
  val hoppingMatrix: PackedSparse[S] = {
    val ret = sparse(2*numLatticeSites, 2*numLatticeSites)
    ret.clear()
    
    val B = 8*math.Pi*B_n / (math.sqrt(3)*w)
    require(w == h) // necessary for B quantization
    
    for (y <- 0 until h;
         x <- 0 until w;
         nn <- 0 until 6;
         (nx, ny) = neighbors(x, y, nn);
         sp <- 0 until 2) {
      val i = matrixIndex(sp, latticeIndex(x, y))
      val j = matrixIndex(sp, latticeIndex(nx, ny))
      
      val (px, py) = position(x, y);
      val (dx, dy) = displacement(x, y, nn);
      val theta = (B/2) * (px*dy - py*dx)
      ret(i, j) = (theta*I).exp * (-t)
    }
    ret.toPacked
  }
}
