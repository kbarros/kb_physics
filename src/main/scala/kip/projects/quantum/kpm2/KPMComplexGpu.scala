package kip.projects.quantum.kpm2

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
import java.util.Arrays
import java.nio.FloatBuffer


class KPMComplexGpu(val cworld: JCudaWorld, H: SparseCsrComplex, s: Int, M: Int, Mq: Int) extends KPMComplex(H, s, M, Mq) {
  val alphaM2 = new Array[Float](2*n*s)
  val alphaM1 = new Array[Float](2*n*s)
  val floatStorage = new ArrayBuf[Float]()
  
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
  
  def convertArrayToSinglePrecision(src: Array[Double], dst: ArrayBuf[Float]) {
    dst.clear()
    for (i <- 0 until src.size) {
      dst.add(src(i).toFloat)
    }
  }
  
  def convertBufferToSinglePrecision(src: ArrayBuf[Double], dst: ArrayBuf[Float]) {
    dst.clear()
    for (i <- 0 until src.size) {
      dst.add(src(i).toFloat)
    }
  }
  
  // Allocate workspace on device appropriate to matrix H
  class State() {
    val indexBytes  = Hs.nnz*Sizeof.INT
    val rowPtrBytes = (n+1)*Sizeof.INT
    val matBytes    = Hs.nnz*2*Sizeof.FLOAT
    val vecBytes    = n*s*2*Sizeof.FLOAT
    
    // Matrix representation on device
    convertBufferToSinglePrecision(Hs.data, floatStorage)
    val cooRowIndex = cworld.allocDeviceArray(Hs.rowIdx.buffer, indexBytes)
    val cooColIndex = cworld.allocDeviceArray(Hs.colIdx.buffer, indexBytes)
    val cooVal      = cworld.allocDeviceArray(floatStorage.buffer, matBytes)
    val csrRowPtr   = cworld.allocDeviceArray(rowPtrBytes)
    // Convert to CSR matrix
    cusparseXcoo2csr(handle, cooRowIndex, Hs.nnz, n, csrRowPtr, CUSPARSE_INDEX_BASE_ZERO)
    
    // Can potentially reorder (dis, djs) indices to improve coalesced memory access
    val dis = Hs.rowIdx.buffer
    val djs = Hs.colIdx.buffer
    val dis_d = cworld.allocDeviceArray(dis, indexBytes)
    val djs_d = cworld.allocDeviceArray(djs, indexBytes)
    
    // Gradient storage, ordered by (dis, djs) indices
    val gradVal_d = cworld.allocDeviceArray(matBytes)
    
    // Random vector storage
    convertArrayToSinglePrecision(R.data, floatStorage)
    val r_d  = cworld.allocDeviceArray(floatStorage.buffer, vecBytes)
    
    // Array storage on device
    val a0_d = cworld.allocDeviceArray(vecBytes)
    val a1_d = cworld.allocDeviceArray(vecBytes)
    val a2_d = cworld.allocDeviceArray(vecBytes)
    val b0_d = cworld.allocDeviceArray(vecBytes)
    val b1_d = cworld.allocDeviceArray(vecBytes)
    val b2_d = cworld.allocDeviceArray(vecBytes)
    
    def deallocate() {
      cworld.freeDeviceArray(cooRowIndex)
      cworld.freeDeviceArray(cooColIndex)
      cworld.freeDeviceArray(cooVal)
      cworld.freeDeviceArray(csrRowPtr)
      cworld.freeDeviceArray(dis_d)
      cworld.freeDeviceArray(djs_d)
      cworld.freeDeviceArray(gradVal_d)
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
      Pointer.to(Array(n)),
      Pointer.to(Array(Hs.nnz)),
      Pointer.to(Array(s)),
      Pointer.to(st.dis_d),
      Pointer.to(st.djs_d),
      Pointer.to(a_d),
      Pointer.to(b_d),
      Pointer.to(Array(scal)),
      Pointer.to(st.gradVal_d));
    val blockSizeX = 64
    val gridSizeX = ((Hs.nnz / blockSizeX) max 1) min 256
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
        n, s, n, // (A.numRows, B.numCols, A.numCols)
        alpha,
        descra, st.cooVal, st.csrRowPtr, st.cooColIndex, // A matrix
        b, n, // (B, B.numRows)
        beta,
        c, n) // (C, C.numRows)
  }
  
  def scaleVector(alpha: Complexd, src_d: Pointer, dst_d: Pointer, st: State) {
    cworld.cpyDeviceToDevice(dst_d, src_d, st.vecBytes)
    cublasCscal(n*s, // number of elements
                cuCmplx(alpha.re.toFloat, alpha.im.toFloat),
                dst_d, 1)
  }
  
  val neg_one = cuCmplx(-1, 0)
  val zero = cuCmplx(0, 0)
  val one  = cuCmplx(1, 0)
  val two  = cuCmplx(2, 0)
  
  def forward(es: EnergyScale) {
    this.es = es
    Hs.fromCsr(H)
    es.scale(Hs)
    val st = new State()
    
    val mu = Array.fill(M)(0.0)
    mu(0) = n           // Tr[T_0[H]] = Tr[1]
    mu(1) = Hs.trace_re // Tr[T_1[H]] = Tr[H]
    
    val a_d = Array(st.a0_d, st.a1_d, st.a2_d)
    
    cworld.cpyDeviceToDevice(a_d(0), st.r_d, st.vecBytes)  // a1 = T_0[H] |r> = 1 |r>
    cgemmH(alpha=one, b=st.r_d, beta=zero, a_d(1), st)     // a2 = T_1[H] |r> = H |r>
    
    for (m <- 2 to M-1) {
      cworld.cpyDeviceToDevice(a_d(2), a_d(0), st.vecBytes)
      cgemmH(alpha=two, b=a_d(1), beta=neg_one, a_d(2), st) // a2 <- alpha_m = T_m[H] r = 2 H a1 - a0 
      
      mu(m) = cublasCdotc(n*s, st.r_d, 1, a_d(2), 1).x      // mu <- r^\dag \dot alpha_2
      
      val temp = a_d(0)
      a_d(0) = a_d(1)
      a_d(1) = a_d(2)
      a_d(2) = temp
    }
    
    gamma = KPMUtil.momentTransform(mu, Mq)
    cworld.cpyDeviceToHost(alphaM2, a_d(0)) // alpha_{M-2}
    cworld.cpyDeviceToHost(alphaM1, a_d(1)) // alpha_{M-1}
    st.deallocate()
  }
  
  def gradient(f: Double=>Double): SparseCsrComplex = {
    val coeff = KPMUtil.expansionCoefficients(M, Mq, f, es)
    val st = new State()
    
    val a_d = Array(st.a0_d, st.a1_d, st.a2_d)
    cworld.cpyHostToDevice(a_d(0), alphaM2) // a0 = alpha_{M-2}
    cworld.cpyHostToDevice(a_d(1), alphaM1) // a1 = alpha_{M-1}
    
    val b_d = Array(st.b0_d, st.b1_d, st.b2_d)
    scaleVector(coeff(M-1), st.r_d, b_d(0), st)  // b0 = c(M-1) r
    cworld.clearDeviceArray(b_d(1), st.vecBytes) // b1 = 0
    
    dX_dH.fromCsr(Hs)
    dX_dH.zero()
    cworld.clearDeviceArray(st.gradVal_d, st.matBytes)
    
    // need special logic since (mu_1) is calculated exactly
    def cp(m: Int): Double = if (m == 1) 0 else coeff(m)
    dX_dH += (coeff(1), 0.0)
    
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
      cublasCaxpy(n*s,                                      // b0 += c(m) r
                  cuCmplx(cp(m).toFloat, 0),
                  st.r_d, 1,
                  b_d(0), 1)
    }
    
    floatStorage.grow(2*dX_dH.nnz)
    cworld.cpyDeviceToHost(floatStorage.buffer, st.gradVal_d, st.matBytes)
    for (idx <- 0 until dX_dH.nnz) {
      val re = floatStorage(2*idx+0)
      val im = floatStorage(2*idx+1)
      dX_dH += (st.dis(idx), st.djs(idx), re, im)
    }
    st.deallocate()
    
    dX_dH *= (1.0/es.mag, 0.0)
    dX_dH
  }
  
  def destroy {
    cublasShutdown()
    // ...
  }
}

