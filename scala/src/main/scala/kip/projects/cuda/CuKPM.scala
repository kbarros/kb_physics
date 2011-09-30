package kip.projects.cuda

import jcuda.cuComplex.cuCmplx
import jcuda.jcusparse.JCusparse.cusparseCcsrmm
import jcuda.jcusparse.JCusparse.cusparseCdotci
import jcuda.jcusparse.JCusparse.cusparseCreate
import jcuda.jcusparse.JCusparse.cusparseCreateMatDescr
import jcuda.jcusparse.JCusparse.cusparseSetMatIndexBase
import jcuda.jcusparse.JCusparse.cusparseSetMatType
import jcuda.jcusparse.JCusparse.cusparseXcoo2csr
import jcuda.jcusparse.cusparseIndexBase.CUSPARSE_INDEX_BASE_ZERO
import jcuda.jcusparse.cusparseMatrixType.CUSPARSE_MATRIX_TYPE_GENERAL
import jcuda.jcusparse.cusparseOperation.CUSPARSE_OPERATION_NON_TRANSPOSE
import jcuda.jcusparse.cusparseHandle
import jcuda.jcusparse.cusparseMatDescr

import kip.projects.quantum._

import smatrix._
import Constructors.complexFlt._
import Scalar.ComplexFlt



object CuKPM extends App {
  val cworld = new JCudaWorld()
  cworld.printDeviceProperties()
  
  val H = sparse(2, 2)
  H(0, 0) = 1
  H(0, 1) = I
  H(1, 0) = 10
  H(1, 1) = 20*I
  val r = fromRows(row(1, 0), row(0, 1))
  val ckpm = new CuKPM(cworld, H.toPacked, null, r, nrand=2, 0)
}


class CuKPM(val cworld: JCudaWorld, val H: PackedSparse[ComplexFlt], val c: Array[Float], val r: Dense[S], val nrand: Int, val seed: Int = 0) {
  val n = H.numRows
   
  // Initialize JCusparse library
  val handle = new cusparseHandle();
  cusparseCreate(handle);
  
  // Create and set up matrix descriptor
  val descra = new cusparseMatDescr();
  cusparseCreateMatDescr(descra);
  cusparseSetMatType(descra, CUSPARSE_MATRIX_TYPE_GENERAL); // CUSPARSE_MATRIX_TYPE_HERMITIAN
  cusparseSetMatIndexBase(descra, CUSPARSE_INDEX_BASE_ZERO);

  val (is, js) = H.definedIndices.unzip
  val nnz = is.size // number of non-zero elements
  println("is " + is)
  println("js " + js)
  // Allocate and copy COO matrix
  val cooRowIndex = cworld.allocDeviceArray(is.toArray)
  val cooColIndex = cworld.allocDeviceArray(js.toArray)
  val cooVal      = cworld.allocDeviceArray(H.data.buffer)
  
  // Convert to CSR matrix
  val csrRowPtr = cworld.allocDeviceArray(new Array[Int](n+1))
  cusparseXcoo2csr(handle, cooRowIndex, nnz, n, csrRowPtr, CUSPARSE_INDEX_BASE_ZERO);
  
  // Dense random vectors
  require(r.numRows == n && r.numCols == nrand)
  val r_d = cworld.allocDeviceArray(r.data.buffer)
  val y_h = Array.fill[Float](n*nrand*2)(0)
  val z_h = Array.fill[Float](n*nrand*2)(0)
  y_h(5) = 1
  val y_d = cworld.allocDeviceArray(y_h)
  val z_d = cworld.allocDeviceArray(z_h)
  
  // Sparse-dense matrix multiply
  cusparseCcsrmm(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                 n, nrand, n, // (A.numRows, B.numCols, A.numCols)
                 cuCmplx(1, 0), // alpha
                 descra, cooVal, csrRowPtr, cooColIndex, // A matrix
                 r_d, n, // (B, B.numRows)
                 cuCmplx(0, 0), // beta
                 z_d, n); // (C, C.numRows)

  // Sparse-dense dot product
  val vecIdxs_d = cworld.allocDeviceArray(Array.tabulate[Int](n*nrand)(identity))
  val dotRes = cuCmplx(0, 0)
  cusparseCdotci(handle,
                 n*nrand, // nnz for dense matrix
                 z_d, vecIdxs_d, // "sparse" vector and (full) indices
                 z_d, // dense vector
                 dotRes, // result
                 CUSPARSE_INDEX_BASE_ZERO // idxBase
                 )
  println("dot product = "+dotRes)
  
  val ret = dense(n, nrand)
  cworld.cpyDeviceToHost(ret.data.buffer, z_d)
  println(ret)
  
  /*
  // Returns: (mu(m), alpha_{M-2}, alpha_{M-1})
  def momentsStochastic(r: Dense[S]): (Array[R], Dense[S], Dense[S]) = {
    val mu = Array.fill(order)(0d)
    mu(0) = n                   // Tr[T_0[H]] = Tr[1]
    mu(1) = H.trace.re          // Tr[T_1[H]] = Tr[H]
    
    val a0 = dense(n, nrand)
    val a1 = dense(n, nrand)
    val a2 = dense(n, nrand)
    
    a0 := r                      // T_0[H] |r> = 1 |r>
    a1 :=* (H, r)                // T_1[H] |r> = H |r>
    
    for (m <- 2 to order-1) {
      a2 := a0; a2.gemm(2, H, a1, -1)  // alpha_m = T_m[H] r = 2 H a1 - a0

      mu(m) = (r dagDot a2).re / nrand
      a0 := a1
      a1 := a2
    }
    
    (mu, a0, a1)
  }
  
  def functionAndGradient(r: Dense[S], c: Array[R], grad: PackedSparse[S]): R = {
    grad.clear()
    
    val a2 = dense(n, nrand)
    val (mu, a0, a1) = momentsStochastic(r)
    
    val b2 = dense(n, nrand)
    val b1 = dense(n, nrand)
    val b0 = r * c(order - 1)
    
    // need special logic since (mu_1) is calculated exactly
    for (i <- 0 until grad.numRows) { grad(i, i) += c(1) }
    def cp(m: Int) = if (m == 1) 0d else c(m)
    
    // cache defined indices for speed
    val (indicesI, indicesJ) = {
      val (i, j) = grad.definedIndices.unzip
      (i.toArray, j.toArray)
    }
    
    for (m <- order-2 to 0 by -1) {
      // a0 = alpha_{m}
      // b0 = beta_{m}

      if (nrand > 1) {
        for ((i, j) <- grad.definedIndices; k <- 0 until nrand) {
          grad(i, j) += (if (m == 0) 1 else 2) * b0(i, k).conj * a0(j, k) / nrand
        }
      }
      // equivalent to above, but much faster. b3 is used as a temporary vector.
      else {
        println("fail; need to cplx conjugate b1 and iterate over nrand")
//        if (m > 1) (b3 :=* (2, b1)) else (b3 := b1)
//        for (iter <- 0 until indicesI.length) {
//          grad.scalar.maddTo(false, b3.data, indicesI(iter), a0.data, indicesJ(iter), grad.data, iter)
//        }
      }
      
      a2 := a1
      b2 := b1
      a1 := a0
      b1 := b0
      a0 := a2; a0.gemm(2, H, a1, -1)                   // a0 = 2 H a1 - a2 
      b0 :=* (cp(m), r); b0.gemm(2, H, b1, 1); b0 -= b2 // b0 = c(m) r + 2 H b1 - b2
    }
    
    (c, mu).zipped.map(_*_).sum
  }
  */
}
