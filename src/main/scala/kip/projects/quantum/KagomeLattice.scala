package kip.projects.quantum

import smatrix._
import ctor._
import kip.math.Vec3


object KagomeLattice extends App {
  import kip.util.Util.{time}
  // time("integrated density")(testIntegratedDensity())
  testEigenvalues()
  
  // Plots the integrated density of states
  def testIntegratedDensity() {
    val q = new KagomeLattice(w=12, h=12, t=1, J_H=0.2, e_min= -10, e_max= 10)
    q.setFieldZeroQOrthogonal(q.field)
    q.fillMatrix(q.matrix)
    
    val H = q.matrix
    require((H - H.dag).norm2.abs < 1e-10, "Found non-hermitian hamiltonian!")
    println("N = "+H.numRows)
    val order = 1000
    val kpm = new KPM(H, nrand=1, seed=0)
    val range = KPM.range(npts=5*order)
    
    val plot = KPM.mkPlot("Integrated density of states")
    KPM.plotLines(plot, (range, KPM.integrateDeltas(range, KPM.eigenvaluesExact(H), moment=0)), "Exact", java.awt.Color.RED)
    KPM.plotLines(plot, (range, KPM.integrate(range, KPM.eigenvaluesApprox(order, range, kpm), moment=0)), "Approx", java.awt.Color.BLACK)
  }
  
  def testEigenvalues() {
    val q = new KagomeLattice(w=16, h=16, t=1, J_H=0.01, e_min= -10, e_max=10)
    val n = q.matrix.numRows
    println("Matrix dim = "+n)
    
    q.setFieldZeroQOrthogonal(q.field)
    q.fillMatrix(q.matrix)
    var eig = KPM.eigenvaluesExact(q.matrix)
    val i_cut = n * 2 / 3
    println(s"Gap between: [${eig(i_cut-2)} ${eig(i_cut-1)}, ${eig(i_cut)} ${eig(i_cut+1)}]")
    println()
    
    def weight(eig: Array[Double], cut: Double) = eig.takeWhile(_ <= cut).map(_.re - cut).sum
    
    val mu = 0.5 * (eig(i_cut-1) + eig(i_cut))
    println(s"Action ${weight(eig, mu)} at mu=$mu")
  }
}


// Notation:
//  d    = vector component (3 dimensional)
//  sp   = Dirac spin index
//  v    = [0,1,2] coordinate of sub-lattice
//  x, y = coordinates on triangular super-lattice
//  nn   = [0,1,2,3] nearest neighbor index (oriented clockwise, starting at 3 o'clock)
//
class KagomeLattice(val w: Int, val h: Int, val t: R, val J_H: R, val e_min: R, val e_max: R) extends KondoHamiltonian {
  require(h % 2 == 0, "Need even number of rows, or hamiltonian is non-hermitian")

  val numLatticeSites = 3*w*h
  
  def coord2idx(v: Int, x: Int, y: Int): Int = {
    v + x*(3) + y*(3*w)
  }
  
  def idx2coord(i: Int) = {
    (i%3, (i/3)%w, (i/(3*w))%h)
  }
  
  def setField(desc: String) {
    desc match {
      case "ferro"    => setFieldFerro(field)
      case "0q_ortho" => setFieldZeroQOrthogonal(field)
    }
  }
  
  def setFieldZeroQOrthogonal(field: Array[R]) {
    for (y <- 0 until h;
         x <- 0 until w;
         v <- 0 until 3) {
      val s = v match {
        case 0 => Seq(1, 0, 0)
        case 1 => Seq(0, 1, 0)
        case 2 => Seq(0, 0, 1)
      }
      for (d <- 0 until 3) { 
        field(fieldIndex(d, coord2idx(v, x, y))) = s(d)
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
  lazy val latticePositions: Array[Vec3] = {
    (for (i <- 0 until numLatticeSites) yield { 
      val (v, x, y) = idx2coord(i)
      val a = 2.0                                 // horizontal distance between letters (A <-> B)
      val b = 0.5*math.sqrt(3.0)*a                // vertical distance between letters   (A <-> D)
      val r = 1 / (2*math.sqrt(3))*a              // distance between letter and number  (A <-> 0)
      val theta = -math.Pi/6 + (2*math.Pi/3)*v    // angle from letter to number
      Vec3(a*x + 0.5*a*y + r*math.cos(theta), b*y + r*math.sin(theta), 0)
    }).toArray
  }
  
  // returns (v, x, y) indices for the `nn`th neighbor site
  def neighbor(v: Int, x: Int, y: Int, nn: Int): (Int, Int, Int) = {
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
      //     - 2 --- 0 -
      case (1, 0) => (2, 0, 1)
      case (1, 1) => (0,-1, 1)
      case (1, 2) => (2, 0, 0)
      case (1, 3) => (0, 0, 0)
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
  
  // in *unscaled* (physical) energy units 
  val hoppingMatrix: PackedSparse[S] = {
    val ret = sparse(2*numLatticeSites, 2*numLatticeSites)
    ret.clear()
        
    for (y <- 0 until h;
         x <- 0 until w;
         v <- 0 until 3;
         nn <- 0 until 4;
         sp <- 0 until 2) {
      val (nv, nx, ny) = neighbor(v, x, y, nn);
      val i = matrixIndex(sp, coord2idx(v, x, y))
      val j = matrixIndex(sp, coord2idx(nv, nx, ny))
      ret(i, j) = -t
    }
    ret.toPacked
  }
}

