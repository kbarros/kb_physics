package kip.projects.cuda

import jcuda.cuComplex
import jcuda.Pointer
import jcuda.Sizeof
import jcuda.cuComplex.cuCmplx
import jcuda.jcublas.JCublas._
import jcuda.jcusparse.JCusparse._
import jcuda.jcusparse.cusparseIndexBase.CUSPARSE_INDEX_BASE_ZERO
import jcuda.jcusparse.cusparseMatrixType.CUSPARSE_MATRIX_TYPE_GENERAL
import jcuda.jcusparse.cusparseOperation.CUSPARSE_OPERATION_NON_TRANSPOSE
import jcuda.jcusparse.cusparseHandle
import jcuda.jcusparse.cusparseMatDescr
import kip.projects.quantum._
import smatrix._
import Constructors.complexFlt._
import Scalar.ComplexFlt
import scala.collection.mutable.ArrayBuffer


object CuKPM extends App {
  val cworld = new JCudaWorld()
  cworld.printDeviceProperties()
  
  val q = new Quantum(w=10, h=10, t=1, J_eff=2, e_min= -10, e_max= 10)  // hopping only: e_min= -6-0.5, e_max= 3+0.5
  val H = q.matrix
  require((H - H.dag).norm2.abs < 1e-10, "Found non-hermitian hamiltonian!")
  println("N = "+H.numRows)

  val order = 5
  val nrand = 20
  val kpm = new KPM(H, order, nrand)
  val r = kpm.randomVector()
  
  
  val ckpm = new CuKPM(cworld, H, order, null, nrand, seed=0)
//  val mu1 = kip.util.Util.time("Cuda")(ckpm.momentsStochastic(r))
//  val (mu2, _, _) = kip.util.Util.time("Cpu")(kpm.momentsStochastic(r))
  
  val c = kpm.expansionCoefficients(de=1e-4, e => e)
  val dH = H.duplicate
  ckpm.functionAndGradient(r, c, dH)
  
//  for (m <- 0 until order) {
//    println("%d = %f : %f".format(m, mu1(m), mu2(m)))
//  }
  
//  val range = kpm.range
//  val plot = KPM.mkPlot("Integrated density of states")
//  KPM.plotLines(plot, (kpm.range, KPM.integrateDeltas(range, KPM.eigenvaluesExact(H), moment=0)), "Exact", java.awt.Color.RED)
//  KPM.plotLines(plot, (kpm.range, KPM.integrate(range, kpm.eigenvaluesApprox(kpm.jacksonKernel), moment=0)), "Approx", java.awt.Color.BLACK)

}


class CuKPM(val cworld: JCudaWorld, val H: PackedSparse[ComplexFlt], val order: Int, val c: Array[Float], val nrand: Int, val seed: Int = 0) {
  val n = H.numRows
  val vecBytes = n*nrand*2*Sizeof.FLOAT
  
  // Initialize JCublas library
  cublasInit()
  
  // Initialize JCusparse library
  val handle = new cusparseHandle();
  cusparseCreate(handle);
  
  // Create and set up matrix descriptor
  val descra = new cusparseMatDescr();
  cusparseCreateMatDescr(descra);
  cusparseSetMatType(descra, CUSPARSE_MATRIX_TYPE_GENERAL); // CUSPARSE_MATRIX_TYPE_HERMITIAN
  cusparseSetMatIndexBase(descra, CUSPARSE_INDEX_BASE_ZERO);

  // Matrix representation on device
  val (is, js) = H.definedIndices.unzip
  val nnz = is.size // number of non-zero elements
  val cooRowIndex = cworld.allocDeviceArray(is.toArray)
  val cooColIndex = cworld.allocDeviceArray(js.toArray)
  val cooVal      = cworld.allocDeviceArray(H.data.buffer)
  val csrRowPtr   = cworld.allocDeviceArray(new Array[Int](n+1))
  // Convert to CSR matrix
  cusparseXcoo2csr(handle, cooRowIndex, nnz, n, csrRowPtr, CUSPARSE_INDEX_BASE_ZERO);
  
  // Diagonal indices, for storing matrix derivative
  val (dis, djs) = {
    val diags_i = Array.fill(n)(new ArrayBuffer[Int]())
    val diags_j = Array.fill(n)(new ArrayBuffer[Int]())
    for ((i, j) <- H.definedIndices) {
      val d = (j - i + n) % n
      diags_i(d) += i
      diags_j(d) += j
    }
//    for (d <- 0 until n;
//         if diags_i(d).size > 0) {
//      println("diag %d, occupancy %d/%d = %g".format(d, diags_i(d).size, n, diags_i(d).size.toDouble / n))
//    }
    (diags_i.flatten, diags_j.flatten)
  }
  require(dis.size == nnz)
  val dis_d = cworld.allocDeviceArray(dis)
  val djs_d = cworld.allocDeviceArray(djs)
  
  // Gradient storage, ordered by diagonal indices
  val gradBytes = nnz*2*Sizeof.FLOAT 
  val gradVal_h = new Array[Float](nnz*2)
  val gradVal_d = cworld.allocDeviceArray(gradVal_h)
  
  // Array of indices for dot product (TODO: delete when cuBlas.cdot available)
  val vecIdxs_d = cworld.allocDeviceArray(Array.tabulate[Int](n*nrand)(identity))
  
  // Array storage on device
  val emptyVec = new Array[Float](n*nrand*2)
  val r_d  = cworld.allocDeviceArray(emptyVec)
  var a0_d = cworld.allocDeviceArray(emptyVec)
  var a1_d = cworld.allocDeviceArray(emptyVec)
  var a2_d = cworld.allocDeviceArray(emptyVec)
  var b0_d = cworld.allocDeviceArray(emptyVec)
  var b1_d = cworld.allocDeviceArray(emptyVec)
  var b2_d = cworld.allocDeviceArray(emptyVec)

  
  def cgemmH(alpha: cuComplex, b: Pointer, beta: cuComplex, c: Pointer) {
    // Sparse-dense matrix multiply
    cusparseCcsrmm(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
        n, nrand, n, // (A.numRows, B.numCols, A.numCols)
        alpha,
        descra, cooVal, csrRowPtr, cooColIndex, // A matrix
        b, n, // (B, B.numRows)
        beta,
        c, n) // (C, C.numRows)
  }
  
  def dotc(x: Pointer, y: Pointer): cuComplex = {
    val ret = cuCmplx(0, 0)
    cusparseCdotci(handle,
        n*nrand, // nnz for dense matrix
        x, vecIdxs_d, // "sparse" vector and (full) indices
        y, // dense vector
        ret, // result
        CUSPARSE_INDEX_BASE_ZERO // idxBase
        )
    ret
  }
  
  def scaleVector(alpha: S#A, src_d: Pointer, dst_d: Pointer) {
    cworld.cpyDeviceToDevice(dst_d, src_d, vecBytes)
    cublasCscal(n*nrand, // number of elemens
                cuCmplx(alpha.re, alpha.im),
                dst_d, 0)
  }
  
  val neg_one = cuCmplx(-1, 0)
  val zero = cuCmplx(0, 0)
  val one  = cuCmplx(1, 0)
  val two  = cuCmplx(2, 0)

  // final two vectors are stored in a0 and a1
  def momentsStochastic(r: Dense[S]): Array[R] = {
    require(r.numRows == n && r.numCols == nrand)

    val mu = Array.fill(order)(0f)
    mu(0) = n                   // Tr[T_0[H]] = Tr[1]
    mu(1) = H.trace.re          // Tr[T_1[H]] = Tr[H]

    cworld.cpyHostToDevice(r_d, r.data.buffer)
    
    cworld.cpyDeviceToDevice(a0_d, r_d, vecBytes) // a1 = T_0[H] |r> = 1 |r>
    cgemmH(alpha=one, b=r_d, beta=zero, a1_d)     // a2 = T_1[H] |r> = H |r>
    
    for (m <- 2 to order-1) {
      cworld.cpyDeviceToDevice(a2_d, a0_d, vecBytes)
      cgemmH(alpha=two, b=a1_d, beta=neg_one, a2_d) // a2 <- alpha_m = T_m[H] r = 2 H a1 - a0 
      
      mu(m) = dotc(r_d, a2_d).x / nrand
      
      val temp = a0_d
      a0_d = a1_d
      a1_d = a2_d
      a2_d = temp
    }
    
    mu
  }

  def functionAndGradient(r: Dense[S], c: Array[R], grad: PackedSparse[S]): R = {
    grad.clear()
    cworld.clearDeviceArray(gradVal_d, gradBytes)

    val mu = momentsStochastic(r) // sets a0_d=alpha_{M-2} and a1_d=alpha_{M-1}
    
    cworld.clearDeviceArray(b1_d, vecBytes) // b1 = 0
    scaleVector(c(order-1), r_d, b0_d)      // b0 = c(order-1) r 

    // need special logic since (mu_1) is calculated exactly
    def cp(m: Int): R = if (m == 1) 0 else c(m)
    for (i <- 0 until grad.numRows) { grad(i, i) += c(1) }
    
    for (m <- order-2 to 0 by -1) {
      // a0 = alpha_{m}
      // b0 = beta_{m}

//      if (nrand > 1) {
//        for ((i, j) <- grad.definedIndices; k <- 0 until nrand) {
//          grad(i, j) += (if (m == 0) 1 else 2) * b0(i, k).conj * a0(j, k) / nrand
//        }
//      }
      
      // (a0, a1, a2) <= (2 H a1 - a2, a0, a1)
      val temp = a2_d
      a2_d = a1_d
      a1_d = a0_d
      a0_d = temp
      cworld.cpyDeviceToDevice(a0_d, a2_d, vecBytes)
      cgemmH(alpha=two, b=a1_d, beta=neg_one, a0_d) // a0 := 2 H a1 - a2 

      // (b0, b1, b2) <= (2 H b1 - b2 + c(m) r, a0, a1)
      val temp2 = b2_d
      b2_d = b1_d
      b1_d = b0_d
      b0_d = temp2
      cworld.cpyDeviceToDevice(b0_d, b2_d, vecBytes)
      cgemmH(alpha=two, b=b1_d, beta=neg_one, b0_d) // b0 := 2 H b1 - b2
      cublasCaxpy(n*nrand,                          // b0 += c(m) r
                  cuCmplx(cp(m), 0),
                  r_d, 0,
                  b0_d, 0)
    }
    
    cworld.cpyDeviceToHost(gradVal_h, gradVal_d)
    for (idx <- 0 until nnz) {
      grad(dis(idx), djs(idx)) += gradVal_h(2*idx+0) + I*gradVal_h(2*idx+1)
    }
    
    val atest = new Array[Float](2*n*nrand)
    cworld.cpyDeviceToHost(atest, a1_d)
    println("PRINTING: ")
    for (i <- 0 until 5) {
      println("%s %s".format(atest(2*i+0) + I*atest(2*i+1), r(i, 0))) 
    }
    
    (c, mu).zipped.map(_*_).sum
  }
  
  def destroy {
    cublasShutdown()
    // ...
  }
}
