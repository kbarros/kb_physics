package kip.projects.quantum

import smatrix._
import Constructors.complexDbl._


object Quantum {
  type T = Scalar.ComplexDbl
}

// Notation:
//  d    = vector component (3 dimensional)
//  sp   = Dirac spin index
//  x, y = coordinates on triangular lattice
//  n    = nearest neighbor index on lattice 
//
class Quantum(w: Int, h: Int, t: Double, J_eff: Double) {
  import Quantum._
  
  val vectorDim = 3
  
  def matrixIndex(sp: Int, x: Int, y: Int): Int = {
    sp + x*(2) + y*(2*w)
  }
  
  def fieldIndex(d: Int, x: Int, y: Int): Int = {
    d + x*(3) + y*(3*w)
  }
  val field: Array[T#A] = Array.fill(vectorDim*w*h)(0) 

  def pauliIndex(sp1: Int, sp2: Int, d: Int): Int = {
    sp1 + sp2*(2) + d*(2*2)
  }
  def pauli = Array[T#A] (
    0, 1,
    1, 0,
    
    0, I, // visually transposed, due to row major ordering
   -I, 0,
    
    1, 0,
    0, -1
  )
  

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
  
  def fillMatrix[M[s <: Scalar] <: Sparse[s, M]](m: M[T]) {
    for (y <- 0 until h;
         x <- 0 until w) {
      
      // neighbor hopping
      val (xn, yn) = neighbors(x, y)
      for (n <- 0 until 6;
           sp <- 0 until 2) {
        val i = matrixIndex(sp, x, y)
        val j = matrixIndex(sp, xn(n), yn(n))
        m(i, j) = -t
      }
      
      // hund coupling 
      for (sp1 <- 0 until 2;
           sp2 <- 0 until 2) {
        
        var coupling = 0: T#A
        for (d <- 0 until 3) {
          coupling += pauli(pauliIndex(d, sp1, sp2)) * field(fieldIndex(d, x, y))
        }
        val i = matrixIndex(sp1, x, y)
        val j = matrixIndex(sp2, x, y)
        m(i, j) = 0 // -J_eff * coupling
      }
    }
  }
  
  
  val A: PackedSparse[Scalar.ComplexDbl] = {
    val ret = sparse(2*h*w, 2*h*w)
    fillMatrix(ret)
    ret.toPacked
  }
}
