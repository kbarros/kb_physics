package kip.math
package linalg

import com.sun.jna.Library
import com.sun.jna.Native
import com.sun.jna.ptr.IntByReference


object Netlib {
  implicit object RealFlt extends NetlibRealFlt[Array[Float]]
  implicit object RealDbl extends NetlibRealDbl[Array[Double]]
  implicit object ComplexFlt extends NetlibComplexFlt[Array[Float]]
  implicit object ComplexDbl extends NetlibComplexDbl[Array[Double]]
  
  lazy val cblas: CblasLib = Native.loadLibrary("vecLib", classOf[CblasLib]).asInstanceOf[CblasLib]
  lazy val lapack: LapackLib = Native.loadLibrary("vecLib", classOf[LapackLib]).asInstanceOf[LapackLib]
}


trait Netlib[@specialized(Float, Double) A, B] {
  // enum CBLAS_ORDER
  val CblasRowMajor=101
  val CblasColMajor=102
    
  // enum CBLAS_TRANSPOSE
  val CblasNoTrans=111
  val CblasTrans=112
  val CblasConjTrans=113
  val AtlasConj=114;

  def gemm(Order: Int, TransA: Int, TransB: Int,
           M: Int, N: Int, K: Int,
           alpha: A,
           A: B, lda: Int,
           B: B, ldb: Int,
           beta: A,
           C: B, ldc: Int)

  def gels(TRANS: String, M: Int, N: Int, NRHS: Int,
           A: B, LDA: Int,
           B: B, LDB: Int,
           WORK: B, LWORK: Int, INFO: IntByReference)

  def geev(JOBVL: String, JOBVR: String,
           N: Int, A: B, LDA: Int,
           WR: B, WI: B,
           VL: B, LDVL: Int,
           VR: B, LDVR: Int,
           WORK: B, LWORK: Int, INFO: IntByReference)

  protected implicit def intToIntByReference(a: Int) = new IntByReference(a)
  protected implicit def complexToDoubleArray(c: Complex) = Array[Double](c.re, c.im)
  protected implicit def complexToFloatArray(c: Complex) = Array[Float](c.re.toFloat, c.im.toFloat)  
}


class NetlibRealFlt[B](implicit s: ScalarData[Float, B] with ScalarBufferedFlt) extends Netlib[Float, B] {
  private implicit def dataToBuffer(data: B) = s.buffer(data)
  import Netlib.cblas._
  import Netlib.lapack._
  
  def gemm(Order: Int, TransA: Int, TransB: Int,
           M: Int, N: Int, K: Int,
           alpha: Float,
           A: B, lda: Int,
           B: B, ldb: Int,
           beta: Float,
           C: B, ldc: Int) =
    cblas_sgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)

  def gels(TRANS: String, M: Int, N: Int, NRHS: Int,
           A: B, LDA: Int,
           B: B, LDB: Int,
           WORK: B, LWORK: Int, INFO: IntByReference) =
    sgels_(TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO)

  def geev(JOBVL: String, JOBVR: String,
           N: Int, A: B, LDA: Int,
           WR: B, WI: B,
           VL: B, LDVL: Int,
           VR: B, LDVR: Int,
           WORK: B, LWORK: Int, INFO: IntByReference) =
    sgeev_(JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
             
}


class NetlibRealDbl[B](implicit s: ScalarData[Double, B] with ScalarBufferedDbl) extends Netlib[Double, B] {
  private implicit def dataToBuffer(data: B) = s.buffer(data)
  import Netlib.cblas._
  import Netlib.lapack._
  
  def gemm(Order: Int, TransA: Int, TransB: Int,
           M: Int, N: Int, K: Int,
           alpha: Double,
           A: B, lda: Int,
           B: B, ldb: Int,
           beta: Double,
           C: B, ldc: Int) =
    cblas_dgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)

  def gels(TRANS: String, M: Int, N: Int, NRHS: Int,
           A: B, LDA: Int,
           B: B, LDB: Int,
           WORK: B, LWORK: Int, INFO: IntByReference) = {
    dgels_(TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO)
  }

  def geev(JOBVL: String, JOBVR: String,
           N: Int, A: B, LDA: Int,
           WR: B, WI: B,
           VL: B, LDVL: Int,
           VR: B, LDVR: Int,
           WORK: B, LWORK: Int, INFO: IntByReference) =
    dgeev_(JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
}


class NetlibComplexFlt[B](implicit s: ScalarData[Complex, B] with ScalarBufferedFlt) extends Netlib[Complex, B] {
  private implicit def dataToBuffer(data: B) = s.buffer(data)
  import Netlib.cblas._
  import Netlib.lapack._
  
  def gemm(Order: Int, TransA: Int, TransB: Int,
           M: Int, N: Int, K: Int,
           alpha: Complex,
           A: B, lda: Int,
           B: B, ldb: Int,
           beta: Complex,
           C: B, ldc: Int) =
    cblas_cgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)

  def gels(TRANS: String, M: Int, N: Int, NRHS: Int,
           A: B, LDA: Int,
           B: B, LDB: Int,
           WORK: B, LWORK: Int, INFO: IntByReference) =
    cgels_(TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO)

  def geev(JOBVL: String, JOBVR: String,
           N: Int, A: B, LDA: Int,
           WR: B, WI: B,
           VL: B, LDVL: Int,
           VR: B, LDVR: Int,
           WORK: B, LWORK: Int, INFO: IntByReference) =
    cgeev_(JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
}


class NetlibComplexDbl[B](implicit s: ScalarData[Complex, B] with ScalarBufferedDbl) extends Netlib[Complex, B] {
  private implicit def dataToBuffer(data: B) = s.buffer(data)
  import Netlib.cblas._
  import Netlib.lapack._
  
  def gemm(Order: Int, TransA: Int, TransB: Int,
           M: Int, N: Int, K: Int,
           alpha: Complex,
           A: B, lda: Int,
           B: B, ldb: Int,
           beta: Complex,
           C: B, ldc: Int) =
    cblas_zgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)

  def gels(TRANS: String, M: Int, N: Int, NRHS: Int,
           A: B, LDA: Int,
           B: B, LDB: Int,
           WORK: B, LWORK: Int, INFO: IntByReference) =
    zgels_(TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO)

  def geev(JOBVL: String, JOBVR: String,
           N: Int, A: B, LDA: Int,
           WR: B, WI: B,
           VL: B, LDVL: Int,
           VR: B, LDVR: Int,
           WORK: B, LWORK: Int, INFO: IntByReference) =
    zgeev_(JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
}



// ------------------------------------------------------------------------------
// BINDINGS
//

trait CblasLib extends Library {
    
  def cblas_sgemm(Order: Int, TransA: Int, TransB: Int,
                  M: Int, N: Int, K: Int,
                  alpha: Float,
                  A: Array[Float], lda: Int,
                  B: Array[Float], ldb: Int,
                  beta: Float,
                  C: Array[Float], ldc: Int)
  def cblas_dgemm(Order: Int, TransA: Int, TransB: Int,
                  M: Int, N: Int, K: Int,
                  alpha: Double,
                  A: Array[Double], lda: Int,
                  B: Array[Double], ldb: Int,
                  beta: Double,
                  C: Array[Double], ldc: Int)
  def cblas_cgemm(Order: Int, TransA: Int, TransB: Int,
                  M: Int, N: Int, K: Int,
                  alpha: Array[Float],
                  A: Array[Float], lda: Int,
                  B: Array[Float], ldb: Int,
                  beta: Array[Float],
                  C: Array[Float], ldc: Int)
  def cblas_zgemm(Order: Int, TransA: Int, TransB: Int,
                  M: Int, N: Int, K: Int,
                  alpha: Array[Double],
                  A: Array[Double], lda: Int,
                  B: Array[Double], ldb: Int,
                  beta: Array[Double],
                  C: Array[Double], ldc: Int)
}


trait LapackLib extends Library {
      
      def sgels_(TRANS: String, M: IntByReference, N: IntByReference, NRHS: IntByReference,
                A: Array[Float], LDA: IntByReference,
                B: Array[Float], LDB: IntByReference,
                WORK: Array[Float], LWORK: IntByReference, INFO: IntByReference)
      def dgels_(TRANS: String, M: IntByReference, N: IntByReference, NRHS: IntByReference,
                A: Array[Double], LDA: IntByReference,
                B: Array[Double], LDB: IntByReference,
                WORK: Array[Double], LWORK: IntByReference, INFO: IntByReference)
      def cgels_(TRANS: String, M: IntByReference, N: IntByReference, NRHS: IntByReference,
                A: Array[Float], LDA: IntByReference,
                B: Array[Float], LDB: IntByReference,
                WORK: Array[Float], LWORK: IntByReference, INFO: IntByReference)
      def zgels_(TRANS: String, M: IntByReference, N: IntByReference, NRHS: IntByReference,
                A: Array[Double], LDA: IntByReference,
                B: Array[Double], LDB: IntByReference,
                WORK: Array[Double], LWORK: IntByReference, INFO: IntByReference)

      def sgeev_(JOBVL: String, JOBVR: String,
                N: IntByReference, A: Array[Float], LDA: IntByReference,
                WR: Array[Float], WI: Array[Float],
                VL: Array[Float], LDVL: IntByReference,
                VR: Array[Float], LDVR: IntByReference,
                WORK: Array[Float], LWORK: IntByReference, INFO: IntByReference)
      def dgeev_(JOBVL: String, JOBVR: String,
                N: IntByReference, A: Array[Double], LDA: IntByReference,
                WR: Array[Double], WI: Array[Double],
                VL: Array[Double], LDVL: IntByReference,
                VR: Array[Double], LDVR: IntByReference,
                WORK: Array[Double], LWORK: IntByReference, INFO: IntByReference)
      def cgeev_(JOBVL: String, JOBVR: String,
                N: IntByReference, A: Array[Float], LDA: IntByReference,
                WR: Array[Float], WI: Array[Float],
                VL: Array[Float], LDVL: IntByReference,
                VR: Array[Float], LDVR: IntByReference,
                WORK: Array[Float], LWORK: IntByReference, INFO: IntByReference)
      def zgeev_(JOBVL: String, JOBVR: String,
                N: IntByReference, A: Array[Double], LDA: IntByReference,
                WR: Array[Double], WI: Array[Double],
                VL: Array[Double], LDVL: IntByReference,
                VR: Array[Double], LDVR: IntByReference,
                WORK: Array[Double], LWORK: IntByReference, INFO: IntByReference)
}
