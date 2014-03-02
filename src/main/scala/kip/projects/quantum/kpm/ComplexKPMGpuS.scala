package kip.projects.quantum.kpm

import jcuda.cuComplex
import jcuda.Pointer
import jcuda.Sizeof
import jcuda.cuComplex.cuCmplx
import jcuda.driver.JCudaDriver._
import jcuda.jcublas.JCublas._
import jcuda.jcusparse.JCusparse._
import jcuda.jcusparse.cusparseIndexBase.CUSPARSE_INDEX_BASE_ZERO
import jcuda.jcusparse.cusparseMatrixType.CUSPARSE_MATRIX_TYPE_GENERAL
import jcuda.jcusparse.cusparseOperation.CUSPARSE_OPERATION_NON_TRANSPOSE
import jcuda.jcusparse.cusparseHandle
import jcuda.jcusparse.cusparseMatDescr
import kip.projects.cuda.CuPreparePtx
import kip.projects.cuda.JCudaWorld
import smatrix._
import smatrix.Constructors.complexDbl._
import scala.collection.mutable.ArrayBuffer
import scala.util.Random


class ComplexKPMGpuS(val cworld: JCudaWorld) extends ComplexKPM {
  import ComplexKPM._
  
  // Initialize JCublas library
  cublasInit()
  
  // Initialize JCusparse library
  val handle = new cusparseHandle()
  cusparseCreate(handle)
  
  // Create and set up matrix descriptor
  val descra = new cusparseMatDescr()
  cusparseCreateMatDescr(descra)
  cusparseSetMatType(descra, CUSPARSE_MATRIX_TYPE_GENERAL)
  cusparseSetMatIndexBase(descra, CUSPARSE_INDEX_BASE_ZERO)
  
  // Allocate workspace on device appropriate to matrix H
  class State(val r: Dense[Cd], val Hs: PackedSparse[Cd]) {
    val n = Hs.numRows
    val s = r.numCols
    val (is, js) = Hs.definedIndices.unzip
    val nnz = is.size // number of non-zero elements
    
    val vecBytes = n*s*2*Sizeof.FLOAT

    // Matrix representation on device
    val cooRowIndex = cworld.allocDeviceArray(is.toArray)
    val cooColIndex = cworld.allocDeviceArray(js.toArray)
    val cooVal      = cworld.allocDeviceArray(Hs.map(_.toComplexf).data.buffer)
    val csrRowPtr   = cworld.allocDeviceArray(new Array[Int](n+1))
    // Convert to CSR matrix
    cusparseXcoo2csr(handle, cooRowIndex, nnz, n, csrRowPtr, CUSPARSE_INDEX_BASE_ZERO)

    // Diagonal indices, for storing matrix derivative
    val (dis, djs) = {
      val diags_i = Array.fill(n)(new ArrayBuffer[Int]())
      val diags_j = Array.fill(n)(new ArrayBuffer[Int]())
      for ((i, j) <- Hs.definedIndices) {
        val d = (j - i + n) % n
            diags_i(d) += i
            diags_j(d) += j
      }
      (diags_i.flatten, diags_j.flatten)
    }
    require(dis.size == nnz)
    val dis_d = cworld.allocDeviceArray(dis)
    val djs_d = cworld.allocDeviceArray(djs)

    // Gradient storage, ordered by diagonal indices
    val gradBytes = nnz*2*Sizeof.FLOAT 
    val gradVal_h = new Array[Float](nnz*2)
    val gradVal_d = cworld.allocDeviceArray(gradVal_h)
    
//    // Array of indices for dot product (TODO: delete when cuBlas.cdot available)
//    val vecIdxs_d = cworld.allocDeviceArray(Array.tabulate[Int](n*s)(identity))
    
    val r_d  = cworld.allocDeviceArray(r.map(_.toComplexf).data.buffer)
    
    // Array storage on device
    val emptyVec = new Array[Float](n*s*2)
    val a0_d = cworld.allocDeviceArray(emptyVec)
    val a1_d = cworld.allocDeviceArray(emptyVec)
    val a2_d = cworld.allocDeviceArray(emptyVec)
    val b0_d = cworld.allocDeviceArray(emptyVec)
    val b1_d = cworld.allocDeviceArray(emptyVec)
    val b2_d = cworld.allocDeviceArray(emptyVec)
    
    def deallocate() {
      cworld.freeDeviceArray(cooRowIndex)
      cworld.freeDeviceArray(cooColIndex)
      cworld.freeDeviceArray(cooVal)
      cworld.freeDeviceArray(csrRowPtr)
      cworld.freeDeviceArray(dis_d)
      cworld.freeDeviceArray(djs_d)
      cworld.freeDeviceArray(gradVal_d)
//      cworld.freeDeviceArray(vecIdxs_d)
      cworld.freeDeviceArray(r_d)
      cworld.freeDeviceArray(a0_d)
      cworld.freeDeviceArray(a1_d)
      cworld.freeDeviceArray(a2_d)
      cworld.freeDeviceArray(b0_d)
      cworld.freeDeviceArray(b1_d)
      cworld.freeDeviceArray(b2_d)
    }
  }
  
  // Kernel function is equivalent to:
  //   for ((i, j) <- grad.definedIndices; k <- 0 until s) {
  //     grad(i, j) += scal * b0(i, k).conj * a0(j, k)
  //   }
  val cudaSource = """
#include <cuComplex.h>
extern "C"
__global__ void accumulateGrad(int n, int nnz, int s, int *dis, int *djs, cuFloatComplex *a, cuFloatComplex *b, float scal, cuFloatComplex *gradVal) {
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  while (idx < nnz) {
    int di = dis[idx];
    int dj = djs[idx];
    float acc_re = 0;
    float acc_im = 0;
    for (int k = 0; k < s; k++) {
      cuFloatComplex cb = b[k*n+di];
      cuFloatComplex ca = a[k*n+dj];
      // acc += conj(b) * a
      acc_re += cuCrealf(cb) * cuCrealf(ca) + cuCimagf(cb) * cuCimagf(ca);
      acc_im += cuCrealf(cb) * cuCimagf(ca) - cuCimagf(cb) * cuCrealf(ca);
    }
    // gradVal[idx] += scal * acc
    gradVal[idx] = cuCaddf(gradVal[idx], make_cuFloatComplex(scal * acc_re, scal * acc_im));
    idx += gridDim.x*blockDim.x;
  }
}"""
  cworld.loadModule(CuPreparePtx.fromSrcString(cudaSource), Seq("accumulateGrad"))
  def accumulateGrad(a_d: Pointer, b_d: Pointer, scal: Float, st: State) {
    val kernelParameters = Pointer.to(
      Pointer.to(Array(st.n)),
      Pointer.to(Array(st.nnz)),
      Pointer.to(Array(st.s)),
      Pointer.to(st.dis_d),
      Pointer.to(st.djs_d),
      Pointer.to(a_d),
      Pointer.to(b_d),
      Pointer.to(Array(scal)),
      Pointer.to(st.gradVal_d));
    val blockSizeX = 64
    val gridSizeX = ((st.nnz / blockSizeX) max 1) min 256
    cuLaunchKernel(cworld.functions("accumulateGrad"),
      gridSizeX, 1, 1, // Grid dimension
      blockSizeX, 1, 1, // Block dimension
      0, null, // Shared memory size and stream
      kernelParameters, null) // Kernel- and extra parameters
    cuCtxSynchronize()
  }
  
  def cgemmH(alpha: cuComplex, b: Pointer, beta: cuComplex, c: Pointer, st: State) {
    // Sparse-dense matrix multiply
    cusparseCcsrmm(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
        st.n, st.s, st.n, // (A.numRows, B.numCols, A.numCols)
        alpha,
        descra, st.cooVal, st.csrRowPtr, st.cooColIndex, // A matrix
        b, st.n, // (B, B.numRows)
        beta,
        c, st.n) // (C, C.numRows)
  }
  
//  def dotc(x: Pointer, y: Pointer, st: State): cuComplex = {
//    val ret = cuCmplx(0, 0)
//    cusparseCdotci(handle,
//        st.n*st.s, // nnz for dense matrix
//        x, st.vecIdxs_d, // "sparse" vector and (full) indices
//        y, // dense vector
//        ret, // result
//        CUSPARSE_INDEX_BASE_ZERO // idxBase
//        )
//    ret
//  }
  
  def scaleVector(alpha: Complexd, src_d: Pointer, dst_d: Pointer, st: State) {
    cworld.cpyDeviceToDevice(dst_d, src_d, st.vecBytes)
    cublasCscal(st.n*st.s, // number of elements
                cuCmplx(alpha.re.toFloat, alpha.im.toFloat),
                dst_d, 1)
  }
  
  val neg_one = cuCmplx(-1, 0)
  val zero = cuCmplx(0, 0)
  val one  = cuCmplx(1, 0)
  val two  = cuCmplx(2, 0)
  
  def forward(M: Int, Mq: Int, r: Dense[Cd], H: PackedSparse[Cd], es: EnergyScale): ForwardData = {
    val st = new State(r, es.scale(H))
    
    val mu = Array.fill(M)(0.0)
    mu(0) = st.n           // Tr[T_0[H]] = Tr[1]
    mu(1) = st.Hs.trace.re // Tr[T_1[H]] = Tr[H]
    
    val a_d = Array(st.a0_d, st.a1_d, st.a2_d)
    
    cworld.cpyDeviceToDevice(a_d(0), st.r_d, st.vecBytes)  // a1 = T_0[H] |r> = 1 |r>
    cgemmH(alpha=one, b=st.r_d, beta=zero, a_d(1), st)     // a2 = T_1[H] |r> = H |r>
    
    for (m <- 2 to M-1) {
      cworld.cpyDeviceToDevice(a_d(2), a_d(0), st.vecBytes)
      cgemmH(alpha=two, b=a_d(1), beta=neg_one, a_d(2), st) // a2 <- alpha_m = T_m[H] r = 2 H a1 - a0 
      
//      mu(m) = dotc(st.r_d, a_d(2), st).x
      mu(m) = cublasCdotc(st.n*st.s, st.r_d, 1, a_d(2), 1).x
      
      val temp = a_d(0)
      a_d(0) = a_d(1)
      a_d(1) = a_d(2)
      a_d(2) = temp
    }
    
    val aM2 = Constructors.complexFlt.dense(st.n, st.s)
    val aM1 = Constructors.complexFlt.dense(st.n, st.s)
    cworld.cpyDeviceToHost(aM2.data.buffer, a_d(0)) // alpha_{M-2}
    cworld.cpyDeviceToHost(aM1.data.buffer, a_d(1)) // alpha_{M-1}
    val gamma = KPMUtil.momentTransform(mu, Mq)
    st.deallocate()
    ForwardData(st.Hs, es, r, mu, gamma, aM2.map(_.toComplexd), aM1.map(_.toComplexd))
  }
  
  def reverse(fd: ForwardData, coeff: Array[Double]): PackedSparse[Cd] = {
    val st = new State(fd.r, fd.Hs)
    val M = coeff.size
    
    val a_d = Array(st.a0_d, st.a1_d, st.a2_d)
    cworld.cpyHostToDevice(a_d(0), fd.aM2.map(_.toComplexf).data.buffer) // a0 = alpha_{M-2}
    cworld.cpyHostToDevice(a_d(1), fd.aM1.map(_.toComplexf).data.buffer) // a1 = alpha_{M-1}
    
    val b_d = Array(st.b0_d, st.b1_d, st.b2_d)
    scaleVector(coeff(M-1), st.r_d, b_d(0), st)  // b0 = c(M-1) r
    cworld.clearDeviceArray(b_d(1), st.vecBytes) // b1 = 0
    
    val dE_dHs = st.Hs.duplicate.clear()
    cworld.clearDeviceArray(st.gradVal_d, st.gradBytes)
    
    // need special logic since (mu_1) is calculated exactly
    def cp(m: Int): Double = if (m == 1) 0 else coeff(m)
    for (i <- 0 until dE_dHs.numRows) { dE_dHs(i, i) += coeff(1) }
    
    for (m <- M-2 to 0 by -1) {
      // a0 = alpha_{m}
      // b0 = beta_{m}
      val scal = (if (m == 0) 1f else 2f)
      accumulateGrad(a_d(0), b_d(0), scal, st)
      
      // (a0, a1, a2) <= (2 H a1 - a2, a0, a1)
      val temp = a_d(2)
      a_d(2) = a_d(1)
      a_d(1) = a_d(0)
      a_d(0) = temp
      cworld.cpyDeviceToDevice(a_d(0), a_d(2), st.vecBytes)
      cgemmH(alpha=two, b=a_d(1), beta=neg_one, a_d(0), st) // a0 := 2 H a1 - a2 
      
      // (b0, b1, b2) <= (2 H b1 - b2 + c(m) r, a0, a1)
      val temp2 = b_d(2)
      b_d(2) = b_d(1)
      b_d(1) = b_d(0)
      b_d(0) = temp2
      cworld.cpyDeviceToDevice(b_d(0), b_d(2), st.vecBytes)
      cgemmH(alpha=two, b=b_d(1), beta=neg_one, b_d(0), st) // b0 := 2 H b1 - b2
      cublasCaxpy(st.n*st.s,                                // b0 += c(m) r
                  cuCmplx(cp(m).toFloat, 0),
                  st.r_d, 1,
                  b_d(0), 1)
    }
    
    cworld.cpyDeviceToHost(st.gradVal_h, st.gradVal_d)
    for (idx <- 0 until st.nnz) {
      dE_dHs(st.dis(idx), st.djs(idx)) += st.gradVal_h(2*idx+0) + I*st.gradVal_h(2*idx+1)
    }
    st.deallocate()
    
    dE_dHs.transform(_ / fd.es.mag)
  }
  
  def destroy {
    cublasShutdown()
    // ...
  }
}
