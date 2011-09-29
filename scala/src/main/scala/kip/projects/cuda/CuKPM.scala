package kip.projects.cuda

import kip.util.Util.{time, notime}
import kip.math.Math._
import kip.enrich._
import smatrix._
import kip.projects.quantum._
import Constructors.complexFlt._
import Scalar.ComplexFlt

import jcuda.jcusparse.JCusparse._
import jcuda.jcusparse.cusparseIndexBase._;
import jcuda.jcusparse.cusparseMatrixType._;
import jcuda.jcusparse.cusparseOperation._
import jcuda.runtime.JCuda._;
import jcuda.runtime.cudaMemcpyKind._;
import jcuda._;
import jcuda.jcusparse._;
import jcuda.runtime.JCuda;
import jcuda.cuComplex._
  
import jcuda.driver.CUdevice_attribute._;
import jcuda.driver.JCudaDriver._;
import java.util._;

import jcuda.driver._;


object CuKPM extends App {
  // Enable exceptions and subsequently omit error checks in this sample
  JCusparse.setExceptionsEnabled(true);
  JCuda.setExceptionsEnabled(true);

  deviceProperties()
  
  // Initialize JCusparse library
  val handle = new cusparseHandle();
  cusparseCreate(handle);

  // Create and set up matrix descriptor
  val descra = new cusparseMatDescr();
  cusparseCreateMatDescr(descra);
  cusparseSetMatType(descra, CUSPARSE_MATRIX_TYPE_GENERAL); // CUSPARSE_MATRIX_TYPE_HERMITIAN
  cusparseSetMatIndexBase(descra, CUSPARSE_INDEX_BASE_ZERO);
  
  val H = sparse(2, 2)
  H(0, 0) = 1
  H(0, 1) = I
  H(1, 0) = 10
  H(1, 1) = 20*I
  val ckpm = new CuKPM(H.toPacked, null, r = col(1, 2), nrand=2, 0)
  
    
  def deviceProperties() {
    cuInit(0);

    // Obtain the number of devices
    val deviceCountArray = Array(0)
    cuDeviceGetCount(deviceCountArray);
    val deviceCount = deviceCountArray(0)
    System.out.println("Found " + deviceCount + " devices");

    for (i <- 0 until deviceCount) {
      val device = new CUdevice();
      cuDeviceGet(device, i);

      // Obtain the device name
      val deviceName = new Array[Byte](1024)
      cuDeviceGetName(deviceName, deviceName.length, device);
      val name = deviceName.map(_.toChar).mkString

      // Obtain the compute capability
      val majorArray = Array(0);
      val minorArray = Array(0);
      cuDeviceComputeCapability(
          majorArray, minorArray, device);
      val major = majorArray(0);
      val minor = minorArray(0);

      println("Device " + i + ": " + name + " with Compute Capability " + major + "." + minor);
    }    
  }
}

class CuKPM(val H: PackedSparse[ComplexFlt], val c: Array[Float], val r: Dense[S], val nrand: Int, val seed: Int = 0) {
  val n = H.numRows
  import CuKPM._
  
  val SizeofComplex = 2 * Sizeof.FLOAT
  val (is, js) = H.definedIndices.unzip
  val nnz = is.size // number of non-zero elements
  println("is " + is)
  println("js " + js)
  // Allocate and copy COO matrix
  val cooRowIndex = new Pointer();
  val cooColIndex = new Pointer();
  val cooVal = new Pointer();
  cudaMalloc(cooRowIndex, nnz * Sizeof.INT);
  cudaMalloc(cooColIndex, nnz * Sizeof.INT);
  cudaMalloc(cooVal,      nnz * SizeofComplex);
  cudaMemcpy(cooRowIndex, Pointer.to(is.toArray),    nnz * Sizeof.INT,    cudaMemcpyHostToDevice);
  cudaMemcpy(cooColIndex, Pointer.to(js.toArray),    nnz * Sizeof.INT,    cudaMemcpyHostToDevice);
  cudaMemcpy(cooVal,      Pointer.to(H.data.buffer), nnz * SizeofComplex, cudaMemcpyHostToDevice);
//  cudaMemcpy(cooVal,      Pointer.to(Array[Float](1, 0, 0, 1)), nnz * SizeofComplex, cudaMemcpyHostToDevice);
  
  // Convert to CSR matrix
  val csrRowPtr = new Pointer();
  cudaMalloc(csrRowPtr, (n+1)*Sizeof.INT);
  cusparseXcoo2csr(handle, cooRowIndex, nnz, n, csrRowPtr, CUSPARSE_INDEX_BASE_ZERO);
  
  // Dense random vectors
  require(r.numRows == n && r.numCols == nrand)
  val r_d = new Pointer();
  cudaMalloc(r_d, n*nrand*SizeofComplex);
  cudaMemcpy(r_d, Pointer.to(r.data.buffer), n*nrand*SizeofComplex, cudaMemcpyHostToDevice)

  val yhost = Array.fill[Float](n*nrand*2)(0)
  yhost(5) = 1
//  yhost(0) = 1
  val y = new Pointer();
  val z = new Pointer();
  cudaMalloc(y, n*nrand*SizeofComplex);
  cudaMalloc(z, n*nrand*SizeofComplex);
  cudaMemcpy(y, Pointer.to(yhost), n*nrand*SizeofComplex, cudaMemcpyHostToDevice)
  cudaMemset(z, 0, n*nrand*SizeofComplex);
  
  // Sparse-dense matrix multiply
  cusparseCcsrmm(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                 n, nrand, n, // (A.numRows, B.numCols, A.numCols)
                 cuCmplx(1, 0), // alpha
                 descra, cooVal, csrRowPtr, cooColIndex, // A matrix
                 r_d, n, // (B, B.numRows)
                 cuCmplx(0, 0), // beta
                 z, n); // (C, C.numRows)

//  val ret = new Array[Float](n*nrand*2)
  val ret = dense(n, nrand)
  cudaMemcpy(Pointer.to(ret.data.buffer), r_d, n*nrand*SizeofComplex, cudaMemcpyDeviceToHost);

//  for (i <- 0 until n*nrand) {
//    println(ret(2*i) + ret(2*i+1)*I)
//  }
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
