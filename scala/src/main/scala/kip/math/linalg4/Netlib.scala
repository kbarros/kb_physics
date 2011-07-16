package kip.math
package linalg4


import java.nio.{FloatBuffer, DoubleBuffer}
import com.sun.jna.Library
import com.sun.jna.Native
import com.sun.jna.ptr.IntByReference

object Netlib {
  implicit object RealFlt    extends NetlibRealFlt
  implicit object RealDbl    extends NetlibRealDbl
  implicit object ComplexFlt extends NetlibComplexFlt
  implicit object ComplexDbl extends NetlibComplexDbl
  
  lazy val cblas: CblasLib = Native.loadLibrary("vecLib", classOf[CblasLib]).asInstanceOf[CblasLib]
  lazy val lapack: LapackLib = Native.loadLibrary("vecLib", classOf[LapackLib]).asInstanceOf[LapackLib]
}


trait Netlib[S <: Scalar] {
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
           alpha: S#A,
           A: S#Buf, lda: Int,
           B: S#Buf, ldb: Int,
           beta: S#A,
           C: S#Buf, ldc: Int)

  def gels(TRANS: String, M: Int, N: Int, NRHS: Int,
           A: S#Buf, LDA: Int,
           B: S#Buf, LDB: Int,
           WORK: S#Buf, LWORK: Int, INFO: IntByReference)

  def geev(JOBVL: String, JOBVR: String,
           N: Int, A: S#Buf, LDA: Int,
           WR: S#Buf, WI: S#Buf,
           VL: S#Buf, LDVL: Int,
           VR: S#Buf, LDVR: Int,
           WORK: S#Buf, LWORK: Int, INFO: IntByReference)

  protected implicit def intToIntByReference(a: Int) = new IntByReference(a)
  protected implicit def complexToFloatBuffer(c: Complex)  = FloatBuffer.wrap (Array[Float](c.re.toFloat, c.im.toFloat))
  protected implicit def complexToDoubleBuffer(c: Complex) = DoubleBuffer.wrap(Array[Double](c.re, c.im))
}


trait NetlibRealFlt extends Netlib[Scalar.RealFlt] {
  import Netlib.cblas._
  import Netlib.lapack._
  
  def gemm(Order: Int, TransA: Int, TransB: Int,
           M: Int, N: Int, K: Int,
           alpha: Float,
           A: FloatBuffer, lda: Int,
           B: FloatBuffer, ldb: Int,
           beta: Float,
           C: FloatBuffer, ldc: Int) =
    cblas_sgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)

  def gels(TRANS: String, M: Int, N: Int, NRHS: Int,
           A: FloatBuffer, LDA: Int,
           B: FloatBuffer, LDB: Int,
           WORK: FloatBuffer, LWORK: Int, INFO: IntByReference) =
    sgels_(TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO)

  def geev(JOBVL: String, JOBVR: String,
           N: Int, A: FloatBuffer, LDA: Int,
           WR: FloatBuffer, WI: FloatBuffer,
           VL: FloatBuffer, LDVL: Int,
           VR: FloatBuffer, LDVR: Int,
           WORK: FloatBuffer, LWORK: Int, INFO: IntByReference) =
    sgeev_(JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
}



class NetlibRealDbl extends Netlib[Scalar.RealDbl] {
  import Netlib.cblas._
  import Netlib.lapack._
  
  def gemm(Order: Int, TransA: Int, TransB: Int,
           M: Int, N: Int, K: Int,
           alpha: Double,
           A: DoubleBuffer, lda: Int,
           B: DoubleBuffer, ldb: Int,
           beta: Double,
           C: DoubleBuffer, ldc: Int) =
    cblas_dgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)

  def gels(TRANS: String, M: Int, N: Int, NRHS: Int,
           A: DoubleBuffer, LDA: Int,
           B: DoubleBuffer, LDB: Int,
           WORK: DoubleBuffer, LWORK: Int, INFO: IntByReference) = {
    dgels_(TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO)
  }

  def geev(JOBVL: String, JOBVR: String,
           N: Int, A: DoubleBuffer, LDA: Int,
           WR: DoubleBuffer, WI: DoubleBuffer,
           VL: DoubleBuffer, LDVL: Int,
           VR: DoubleBuffer, LDVR: Int,
           WORK: DoubleBuffer, LWORK: Int, INFO: IntByReference) =
    dgeev_(JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
}

class NetlibComplexFlt extends Netlib[Scalar.ComplexFlt] {
  import Netlib.cblas._
  import Netlib.lapack._
  
  def gemm(Order: Int, TransA: Int, TransB: Int,
           M: Int, N: Int, K: Int,
           alpha: Complex,
           A: FloatBuffer, lda: Int,
           B: FloatBuffer, ldb: Int,
           beta: Complex,
           C: FloatBuffer, ldc: Int) =
    cblas_cgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)

  def gels(TRANS: String, M: Int, N: Int, NRHS: Int,
           A: FloatBuffer, LDA: Int,
           B: FloatBuffer, LDB: Int,
           WORK: FloatBuffer, LWORK: Int, INFO: IntByReference) =
    cgels_(TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO)

  def geev(JOBVL: String, JOBVR: String,
           N: Int, A: FloatBuffer, LDA: Int,
           WR: FloatBuffer, WI: FloatBuffer,
           VL: FloatBuffer, LDVL: Int,
           VR: FloatBuffer, LDVR: Int,
           WORK: FloatBuffer, LWORK: Int, INFO: IntByReference) =
    cgeev_(JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
}


class NetlibComplexDbl extends Netlib[Scalar.ComplexDbl] {
  import Netlib.cblas._
  import Netlib.lapack._
  
  def gemm(Order: Int, TransA: Int, TransB: Int,
           M: Int, N: Int, K: Int,
           alpha: Complex,
           A: DoubleBuffer, lda: Int,
           B: DoubleBuffer, ldb: Int,
           beta: Complex,
           C: DoubleBuffer, ldc: Int) =
    cblas_zgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)

  def gels(TRANS: String, M: Int, N: Int, NRHS: Int,
           A: DoubleBuffer, LDA: Int,
           B: DoubleBuffer, LDB: Int,
           WORK: DoubleBuffer, LWORK: Int, INFO: IntByReference) =
    zgels_(TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO)

  def geev(JOBVL: String, JOBVR: String,
           N: Int, A: DoubleBuffer, LDA: Int,
           WR: DoubleBuffer, WI: DoubleBuffer,
           VL: DoubleBuffer, LDVL: Int,
           VR: DoubleBuffer, LDVR: Int,
           WORK: DoubleBuffer, LWORK: Int, INFO: IntByReference) =
    zgeev_(JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
}



// ------------------------------------------------------------------------------
// BINDINGS
//

trait CblasLib extends Library {
    
  def cblas_sgemm(Order: Int, TransA: Int, TransB: Int,
                  M: Int, N: Int, K: Int,
                  alpha: Float,
                  A: FloatBuffer, lda: Int,
                  B: FloatBuffer, ldb: Int,
                  beta: Float,
                  C: FloatBuffer, ldc: Int)
  def cblas_dgemm(Order: Int, TransA: Int, TransB: Int,
                  M: Int, N: Int, K: Int,
                  alpha: Double,
                  A: DoubleBuffer, lda: Int,
                  B: DoubleBuffer, ldb: Int,
                  beta: Double,
                  C: DoubleBuffer, ldc: Int)
  def cblas_cgemm(Order: Int, TransA: Int, TransB: Int,
                  M: Int, N: Int, K: Int,
                  alpha: FloatBuffer,
                  A: FloatBuffer, lda: Int,
                  B: FloatBuffer, ldb: Int,
                  beta: FloatBuffer,
                  C: FloatBuffer, ldc: Int)
  def cblas_zgemm(Order: Int, TransA: Int, TransB: Int,
                  M: Int, N: Int, K: Int,
                  alpha: DoubleBuffer,
                  A: DoubleBuffer, lda: Int,
                  B: DoubleBuffer, ldb: Int,
                  beta: DoubleBuffer,
                  C: DoubleBuffer, ldc: Int)
}


trait LapackLib extends Library {
      
      def sgels_(TRANS: String, M: IntByReference, N: IntByReference, NRHS: IntByReference,
                A: FloatBuffer, LDA: IntByReference,
                B: FloatBuffer, LDB: IntByReference,
                WORK: FloatBuffer, LWORK: IntByReference, INFO: IntByReference)
      def dgels_(TRANS: String, M: IntByReference, N: IntByReference, NRHS: IntByReference,
                A: DoubleBuffer, LDA: IntByReference,
                B: DoubleBuffer, LDB: IntByReference,
                WORK: DoubleBuffer, LWORK: IntByReference, INFO: IntByReference)
      def cgels_(TRANS: String, M: IntByReference, N: IntByReference, NRHS: IntByReference,
                A: FloatBuffer, LDA: IntByReference,
                B: FloatBuffer, LDB: IntByReference,
                WORK: FloatBuffer, LWORK: IntByReference, INFO: IntByReference)
      def zgels_(TRANS: String, M: IntByReference, N: IntByReference, NRHS: IntByReference,
                A: DoubleBuffer, LDA: IntByReference,
                B: DoubleBuffer, LDB: IntByReference,
                WORK: DoubleBuffer, LWORK: IntByReference, INFO: IntByReference)

      def sgeev_(JOBVL: String, JOBVR: String,
                N: IntByReference, A: FloatBuffer, LDA: IntByReference,
                WR: FloatBuffer, WI: FloatBuffer,
                VL: FloatBuffer, LDVL: IntByReference,
                VR: FloatBuffer, LDVR: IntByReference,
                WORK: FloatBuffer, LWORK: IntByReference, INFO: IntByReference)
      def dgeev_(JOBVL: String, JOBVR: String,
                N: IntByReference, A: DoubleBuffer, LDA: IntByReference,
                WR: DoubleBuffer, WI: DoubleBuffer,
                VL: DoubleBuffer, LDVL: IntByReference,
                VR: DoubleBuffer, LDVR: IntByReference,
                WORK: DoubleBuffer, LWORK: IntByReference, INFO: IntByReference)
      def cgeev_(JOBVL: String, JOBVR: String,
                N: IntByReference, A: FloatBuffer, LDA: IntByReference,
                WR: FloatBuffer, WI: FloatBuffer,
                VL: FloatBuffer, LDVL: IntByReference,
                VR: FloatBuffer, LDVR: IntByReference,
                WORK: FloatBuffer, LWORK: IntByReference, INFO: IntByReference)
      def zgeev_(JOBVL: String, JOBVR: String,
                N: IntByReference, A: DoubleBuffer, LDA: IntByReference,
                WR: DoubleBuffer, WI: DoubleBuffer,
                VL: DoubleBuffer, LDVL: IntByReference,
                VR: DoubleBuffer, LDVR: IntByReference,
                WORK: DoubleBuffer, LWORK: IntByReference, INFO: IntByReference)
}

