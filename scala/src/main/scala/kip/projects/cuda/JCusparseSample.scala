package kip.projects.cuda

/*
 * JCusparse - Java bindings for CUSPARSE, the NVIDIA CUDA sparse
 * matrix library, to be used with JCuda
 *
 * Copyright (c) 2010 Marco Hutter - http://www.jcuda.org
 */

import jcuda.jcusparse.JCusparse._
import jcuda.jcusparse.cusparseIndexBase.CUSPARSE_INDEX_BASE_ZERO;
import jcuda.jcusparse.cusparseMatrixType.CUSPARSE_MATRIX_TYPE_GENERAL;
import jcuda.jcusparse.cusparseOperation.CUSPARSE_OPERATION_NON_TRANSPOSE;
import jcuda.runtime.JCuda._;
import jcuda.runtime.cudaMemcpyKind._;
import jcuda._;
import jcuda.jcusparse._;
import jcuda.runtime.JCuda;

/**
 * A sample application showing how to use JCusparse.<br />
 * <br />
 * This sample has been ported from the NVIDIA CURAND
 * documentation example.
 */

object JCusparseSample extends App {
  // Enable exceptions and subsequently omit error checks in this sample
  JCusparse.setExceptionsEnabled(true);
  JCuda.setExceptionsEnabled(true);

  System.out.println("Testing example");

  // Create the following sparse test matrix in COO format
  // | 1.0     2.0 3.0 |
  // |     4.0 2.0 3.0 |
  // | 5.0     6.0 7.0 |
  // |     8.0     9.0 |
  val n = 4;
  val nnz = 9;
  val cooRowIndexHostPtr = new Array[Int](nnz)
  val cooColIndexHostPtr = new Array[Int](nnz)
  val cooValHostPtr = new Array[Float](nnz)

  cooRowIndexHostPtr(0) = 0; cooColIndexHostPtr(0) = 0; cooValHostPtr(0) = 1.0f;
  cooRowIndexHostPtr(1) = 0; cooColIndexHostPtr(1) = 2; cooValHostPtr(1) = 2.0f;
  cooRowIndexHostPtr(2) = 0; cooColIndexHostPtr(2) = 3; cooValHostPtr(2) = 3.0f;
  cooRowIndexHostPtr(3) = 1; cooColIndexHostPtr(3) = 1; cooValHostPtr(3) = 4.0f;
  cooRowIndexHostPtr(4) = 2; cooColIndexHostPtr(4) = 0; cooValHostPtr(4) = 5.0f;
  cooRowIndexHostPtr(5) = 2; cooColIndexHostPtr(5) = 2; cooValHostPtr(5) = 6.0f;
  cooRowIndexHostPtr(6) = 2; cooColIndexHostPtr(6) = 3; cooValHostPtr(6) = 7.0f;
  cooRowIndexHostPtr(7) = 3; cooColIndexHostPtr(7) = 1; cooValHostPtr(7) = 8.0f;
  cooRowIndexHostPtr(8) = 3; cooColIndexHostPtr(8) = 3; cooValHostPtr(8) = 9.0f;

  // Print the matrix
  println("Input data:\n");
  for (i <- 0 until nnz) {
    println("cooRowIndedHostPtr(%d)=%d  ", i, cooRowIndexHostPtr(i));
    println("cooColIndedHostPtr(%d)=%d  ", i, cooColIndexHostPtr(i));
    println("cooValHostPtr(%d)=%f     \n", i, cooValHostPtr(i));
  }

  // Create a sparse and a dense vector
  // xVal=(100.0, 200.0, 400.0) (sparse)
  // xInd=(0      1      3    )
  // y   =(10.0, 20.0, 30.0, 40.0 | 50.0, 60.0, 70.0, 80.0) (dense)
  val nnz_vector = 3;
  val xIndHostPtr = new Array[Int](nnz_vector);
  val xValHostPtr = new Array[Float](nnz_vector);
  val yHostPtr = new Array[Float](2 * n);
  val zHostPtr = new Array[Float](2 * (n + 1));

  yHostPtr(0) = 10.0f; xIndHostPtr(0) = 0; xValHostPtr(0) = 100.0f;
  yHostPtr(1) = 20.0f; xIndHostPtr(1) = 1; xValHostPtr(1) = 200.0f;
  yHostPtr(2) = 30.0f;
  yHostPtr(3) = 40.0f; xIndHostPtr(2) = 3; xValHostPtr(2) = 400.0f;
  yHostPtr(4) = 50.0f;
  yHostPtr(5) = 60.0f;
  yHostPtr(6) = 70.0f;
  yHostPtr(7) = 80.0f;

  // Print the vectors
  for (j <- 0 until 2) {
    for (i <- 0 until n) {
      println("yHostPtr(%d,%d)=%f\n", i, j, yHostPtr(i + n * j));
    }
  }
  for (i <- 0 until nnz_vector) {
    println("xIndHostPtr(%d)=%d  ", i, xIndHostPtr(i));
    println("xValHostPtr(%d)=%f\n", i, xValHostPtr(i));
  }

  // Allocate GPU memory and copy the matrix and vectors into it
  val cooRowIndex = new Pointer();
  val cooColIndex = new Pointer();
  val cooVal = new Pointer();
  val xInd = new Pointer();
  val xVal = new Pointer();
  val y = new Pointer();

  cudaMalloc(cooRowIndex, nnz * Sizeof.INT);
  cudaMalloc(cooColIndex, nnz * Sizeof.INT);
  cudaMalloc(cooVal, nnz * Sizeof.FLOAT);
  cudaMalloc(y, 2 * n * Sizeof.FLOAT);
  cudaMalloc(xInd, nnz_vector * Sizeof.INT);
  cudaMalloc(xVal, nnz_vector * Sizeof.FLOAT);
  cudaMemcpy(cooRowIndex, Pointer.to(cooRowIndexHostPtr), nnz * Sizeof.INT, cudaMemcpyHostToDevice);
  cudaMemcpy(cooColIndex, Pointer.to(cooColIndexHostPtr), nnz * Sizeof.INT, cudaMemcpyHostToDevice);
  cudaMemcpy(cooVal, Pointer.to(cooValHostPtr), nnz * Sizeof.FLOAT, cudaMemcpyHostToDevice);
  cudaMemcpy(y, Pointer.to(yHostPtr), 2 * n * Sizeof.FLOAT, cudaMemcpyHostToDevice);
  cudaMemcpy(xInd, Pointer.to(xIndHostPtr), nnz_vector * Sizeof.INT, cudaMemcpyHostToDevice);
  cudaMemcpy(xVal, Pointer.to(xValHostPtr), nnz_vector * Sizeof.FLOAT, cudaMemcpyHostToDevice);

  // Initialize JCusparse library
  val handle = new cusparseHandle();
  cusparseCreate(handle);

  // Create and set up matrix descriptor
  val descra = new cusparseMatDescr();
  cusparseCreateMatDescr(descra);
  cusparseSetMatType(descra, CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparseSetMatIndexBase(descra, CUSPARSE_INDEX_BASE_ZERO);

  // Exercise conversion routines (convert matrix from COO 2 CSR format)
  val csrRowPtr = new Pointer();
  cudaMalloc(csrRowPtr, (n + 1) * Sizeof.INT);
  cusparseXcoo2csr(handle, cooRowIndex, nnz, n,
    csrRowPtr, CUSPARSE_INDEX_BASE_ZERO);
  //csrRowPtr = (0 3 4 7 9)

  // Exercise Level 1 routines (scatter vector elements)
  val yn = y.withByteOffset(n * Sizeof.FLOAT);
  cusparseSsctr(handle, nnz_vector, xVal, xInd,
    yn, CUSPARSE_INDEX_BASE_ZERO);
  // y = (10 20 30 40 | 100 200 70 400)

  // Exercise Level 2 routines (csrmv)
  val y0 = y.withByteOffset(0);
  cusparseScsrmv(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, n, n, 2.0f,
    descra, cooVal, csrRowPtr, cooColIndex, y0, 3.0f, yn);

  // Print intermediate results (y)
  // y = (10 20 30 40 | 680 760 1230 2240)
  cudaMemcpy(Pointer.to(yHostPtr), y, 2 * n * Sizeof.FLOAT, cudaMemcpyDeviceToHost);
  for (j <- 0 until 2) {
    for (i <- 0 until n) {
      println("yHostPtr(%d,%d)=%f\n", i, j, yHostPtr(i + n * j));
    }
  }

  // Exercise Level 3 routines (csrmm)
  val z = new Pointer();
  cudaMalloc(z, 2 * (n + 1) * Sizeof.FLOAT);
  cudaMemset(z, 0, 2 * (n + 1) * Sizeof.FLOAT);
  cusparseScsrmm(handle, CUSPARSE_OPERATION_NON_TRANSPOSE, n, 2, n,
    5.0f, descra, cooVal, csrRowPtr, cooColIndex, y, n, 0.0f, z, n + 1);

  // Print final results (z)
  // z = (950 400 2550 2600 0 | 49300 15200 132300 131200 0)
  cudaMemcpy(Pointer.to(zHostPtr), z, 2 * (n + 1) * Sizeof.FLOAT, cudaMemcpyDeviceToHost);

  System.out.printf("Final results:\n");
  for (j <- 0 until 2) {
    for (i <- 0 until n - 1) {
      println("z(%d,%d)=%f\n", i, j, zHostPtr(i + (n + 1) * j));
    }
  }
  if ((zHostPtr(0) != 950.0) ||
    (zHostPtr(1) != 400.0) ||
    (zHostPtr(2) != 2550.0) ||
    (zHostPtr(3) != 2600.0) ||
    (zHostPtr(4) != 0.0) ||
    (zHostPtr(5) != 49300.0) ||
    (zHostPtr(6) != 15200.0) ||
    (zHostPtr(7) != 132300.0) ||
    (zHostPtr(8) != 131200.0) ||
    (zHostPtr(9) != 0.0) ||
    (yHostPtr(0) != 10.0) ||
    (yHostPtr(1) != 20.0) ||
    (yHostPtr(2) != 30.0) ||
    (yHostPtr(3) != 40.0) ||
    (yHostPtr(4) != 680.0) ||
    (yHostPtr(5) != 760.0) ||
    (yHostPtr(6) != 1230.0) ||
    (yHostPtr(7) != 2240.0)) {
    System.out.println("example test FAILED");
  } else {
    System.out.println("example test PASSED");
  }

  // Clean up
  cudaFree(y);
  cudaFree(z);
  cudaFree(xInd);
  cudaFree(xVal);
  cudaFree(csrRowPtr);
  cudaFree(cooRowIndex);
  cudaFree(cooColIndex);
  cudaFree(cooVal);
  cusparseDestroy(handle);
}
