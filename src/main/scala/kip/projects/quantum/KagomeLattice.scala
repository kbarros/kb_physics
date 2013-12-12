package kip.projects.quantum

import smatrix._
import ctor._


object KagomeLattice extends App {
  testIntegratedDensity()
  
//  time("eigenvalues")(testEigenvalues())
  
  // Plots the integrated density of states
  def testIntegratedDensity() {
    val q = new KagomeLattice(w=16, h=16, t=1, J_H=0, B_n=0, e_min= -10, e_max= 10)
    q.setFieldQZeroOrthogonal(q.field)
    q.fillMatrix(q.matrix)
    
    val H = q.matrix
    require((H - H.dag).norm2.abs < 1e-10, "Found non-hermitian hamiltonian!")
    println("N = "+H.numRows)
    val order = 500
    val kpm = new KPM(H, nrand=1, seed=0)
    val range = KPM.range(npts=5*order)
    
    val plot = KPM.mkPlot("Integrated density of states")
    KPM.plotLines(plot, (range, KPM.integrateDeltas(range, KPM.eigenvaluesExact(H), moment=0)), "Exact", java.awt.Color.RED)
    KPM.plotLines(plot, (range, KPM.integrate(range, KPM.eigenvaluesApprox(order, range, kpm), moment=0)), "Approx", java.awt.Color.BLACK)
  }
}


//class KondoHamiltonian {
//}


// Notation:
//  d    = vector component (3 dimensional)
//  sp   = Dirac spin index
//  v    = [0,1,2] coordinate on sublattice
//  x, y = coordinates on triangular lattice
//
//  nn   = [0,1,2,3] nearest neighbor index on lattice 
//
class KagomeLattice(val w: Int, val h: Int, val t: R, val J_H: R, val B_n: Int, val e_min: R, val e_max: R) {
  require(h % 2 == 0, "Need even number of rows, or hamiltonian is non-hermitian")
  val vectorDim = 3
  
  val e_avg   = (e_max + e_min)/2
  val e_scale = (e_max - e_min)/2
  
  def matrixIndex(sp: Int, v: Int, x: Int, y: Int): Int = {
    sp + v*(2) + x*(2*3) + y*(2*3*w)
  }
  
  def fieldIndex(d: Int, v: Int, x: Int, y: Int): Int = {
    d + v*(3) + x*(3*3) + y*(3*3*w)
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
    val ret = new Array[R](3*3*w*h)
    setFieldFerro(ret)
    ret
  }
  val delField: Array[R] = Array.fill(3*3*w*h)(0)
  
  val hoppingMatrix: PackedSparse[S] = {
    val ret = sparse(2*3*w*h, 2*3*w*h)
    fillHoppingMatrix(ret)
    ret.toPacked
  }

  val matrix = {
    val ret = sparse(2*3*w*h, 2*3*w*h): HashSparse[S]
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
    field.transform(_ => rand.nextGaussian())
    normalizeField(field)
  }
  
  def setFieldQZeroOrthogonal(field: Array[R]) {
    for (y <- 0 until h;
         x <- 0 until w;
         v <- 0 until 3) {
      val s = v match {
        case 0 => Seq(1, 0, 0)
        case 1 => Seq(0, 1, 0)
        case 2 => Seq(0, 0, 1)
      }
      for (d <- 0 until vectorDim) { 
        field(fieldIndex(d, v, x, y)) = s(d)
      }
    }
    normalizeField(field)
  }
  
  def normalizeField(field: Array[R], validate: Boolean = false) {
    for (y <- 0 until h;
         x <- 0 until w;
         v <- 0 until 3) {
      var acc = 0d
      for (d <- 0 until 3) {
        acc += field(fieldIndex(d, v, x, y)).abs2
      }
      acc = math.sqrt(acc)
      if (validate && !(acc > 0.95 && acc < 1.05))
        println("Vector magnitude %g deviates too far from normalization".format(acc))
      for (d <- 0 until 3) {
        field(fieldIndex(d, v, x, y)) /= acc
      }
    }
  }
  
  // remove component of dS that is parallel to field S
  def projectTangentField(S: Array[R], dS: Array[R]) {
    for (y <- 0 until h;
         x <- 0 until w;
         v <- 0 until 3) {
      var s_dot_s = 0d
      var s_dot_ds = 0d
      for (d <- 0 until 3) {
        val i = fieldIndex(d, v, x, y)
        s_dot_s  += S(i)*S(i)
        s_dot_ds += S(i)*dS(i)
      }
      val alpha = s_dot_ds / s_dot_s
      for (d <- 0 until 3) {
        val i = fieldIndex(d, v, x, y)
        dS(i) -= alpha * S(i)
      }
    }
  }

  //
  //         1         1         1
  //        /D\       /E\       /F\
  //   --- 2 - 0 --- 2 - 0 --- 2 - 0
  //   \ /       \ /       \ /
  //    1         1         1
  //   /A\       /B\       /C\
  //  2 - 0 --- 2 - 0 --- 2 - 0 ---
  //        \ /       \ /       \ /
  //
  def position(v: Int, x: Int, y: Int): (Double, Double) = {
    val a = 1                                   // horizontal distance between letters (A <-> B)
    val b = 0.5*math.sqrt(3.0)*a                // vertical distance between letters   (A <-> D)
	val r = 1 / (2*math.sqrt(3))*a              // distance between letter and number  (A <-> 0)
	val theta = -math.Pi/6 + (2*math.Pi/3)*v    // angle from letter to number
	(a*x + 0.5*a*y + r*math.cos(theta), b*y + r*math.sin(theta))
  }
    
  // returns (v, x, y) indices for the `nn`th neighbor site
  def neighbors(v: Int, x: Int, y: Int, nn: Int): (Int, Int, Int) = {
	val (nv, xdel, ydel) = (v, nn) match {

    //     1        
    //    /D\       /E
	//   2 - 0 --- 2 - 
	//         \ /
	//          1
	//         /B\
	case (0, 0) => (2, 1, 0)
	case (0, 1) => (1, 0, 0)
	case (0, 2) => (2, 0, 0)
	case (0, 3) => (1, 1,-1)

    //     D\       /E
	//     - 0 --- 2 - 
	//         \ /  
	//          1   
	//         /B\  
	//     - 0 --- 2 -
	case (1, 0) => (2, 0, 1)
	case (1, 1) => (0,-1, 1)
	case (1, 2) => (0, 0, 0)
	case (1, 3) => (2, 0, 0)

	//              1   
	//    D\       /E\  
	//    - 0 --- 2 - 0 -
	//        \ /       
	//         1        
	//        /B\
	case (2, 0) => (0, 0, 0)
	case (2, 1) => (1, 0, 0)
	case (2, 2) => (0,-1, 0)
	case (2, 3) => (1, 0,-1)
	}
	(nv, (x+xdel+w)%w, (y+ydel+h)%h)
  }
    
  // returns (dx, dy) positions
  def displacement(v: Int, x: Int, y: Int, nn: Int): (Double, Double) = {
    val (x0, y0) = position(v, x, y)
    val (nv, nx, ny) = neighbors(v, x, y, nn)
    val (x1, y1) = position(nv, nx, ny)
    (x1-x0, y1-y0)
  }
  
  def fillHoppingMatrix[M[s <: Scalar] <: Sparse[s, M]](m: M[S]) {
    m.clear()
    
    // val B = 8*math.Pi*B_n / (math.sqrt(3)*w)
    require(w == h) // necessary for B quantization
    
    for (y <- 0 until h;
         x <- 0 until w;
         v <- 0 until 3;
         nn <- 0 until 4;
         sp <- 0 until 2) {
      val (nv, nx, ny) = neighbors(v, x, y, nn);
      val i = matrixIndex(sp, v, x, y)
      val j = matrixIndex(sp, nv, nx, ny)
      
      // val (px, py) = position(v, x, y);
      // val (dx, dy) = displacement(v, x, y, nn);
      // val theta = (B/2) * (px*dy - py*dx)
      // m(i, j) = (theta*I).exp * (-t)
      m(i, j) = -t
    }
  }
  
  def fillMatrix[M[s <: Scalar] <: Sparse[s, M]](m: M[S]) {
    m.clear()
    
    for ((i, j) <- hoppingMatrix.definedIndices) {
      m(i, j) += hoppingMatrix(i, j)
    }
    
    // loop over all lattice sites
    for (y <- 0 until h;
         x <- 0 until w;
         v <- 0 until 3) {
      
      // hund coupling 
      for (sp1 <- 0 until 2;
           sp2 <- 0 until 2) {
        
        var coupling = 0: S#A
        for (d <- 0 until 3) {
          coupling += pauli(pauliIndex(sp1, sp2, d)) * field(fieldIndex(d, v, x, y))
        }
        val i = matrixIndex(sp1, v, x, y)
        val j = matrixIndex(sp2, v, x, y)
        m(i, j) = -J_H * coupling
      }
    }
    
    // scale matrix appropriately so that eigenvalues lie between -1 and +1
    for (i <- 0 until m.numRows) { m(i,i) -= e_avg }
    m /= e_scale

//    // Make sure hamiltonian is hermitian
//    val H = m.toDense
//    require((H - H.dag).norm2.abs < 1e-6, "Found non-hermitian hamiltonian!")
  }

  // convert from physical energy units to scaled energy units
  def scaleEnergy(x: R): R = {
    (x - e_avg) / e_scale
  }
  
  // Use chain rule to transform derivative wrt matrix elements dF/dH, into derivative wrt spin indices
  //   dF/dS = dF/dH dH/dS
  // In both factors, H is assumed to be dimensionless (scaled energy). If F is also dimensionless, then it
  // may be desired to multiply the final result by the energy scale.
  def fieldDerivative(dFdH: PackedSparse[S], dFdS: Array[R]) {
    // loop over all lattice sites and vector indices
    for (y <- 0 until h;
         x <- 0 until w;
         v <- 0 until 3;
         d <- 0 until 3) {
      
      var dCoupling: S#A = 0
      for (sp1 <- 0 until 2;
           sp2 <- 0 until 2) {
        val i = matrixIndex(sp1, v, x, y)
        val j = matrixIndex(sp2, v, x, y)
        dCoupling += dFdH(i, j) * pauli(pauliIndex(sp1, sp2, d))
      }
      require(math.abs(dCoupling.im) < 1e-5, "Imaginary part of field derivative non-zero: " + dCoupling.im)
      dFdS(fieldIndex(d, v, x, y)) = -J_H * dCoupling.re
    }
    
    // the derivative is perpendicular to the direction of S, due to constraint |S|=1 
    projectTangentField(field, dFdS)
    
    // properly scale the factor dH/dS, corresponding to scaled H
    dFdS.transform(_ / e_scale)
  }
}
