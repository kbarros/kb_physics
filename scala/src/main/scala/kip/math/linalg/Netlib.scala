package kip.math.linalg

import com.sun.jna.Library
import com.sun.jna.ptr.IntByReference


object CblasLib {
  // enum CBLAS_ORDER
  val CblasRowMajor=101
  val CblasColMajor=102
    
  // enum CBLAS_TRANSPOSE
  val CblasNoTrans=111
  val CblasTrans=112
  val CblasConjTrans=113
  val AtlasConj=114;
}

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
      
      def sgels_(TRANS: String, M: Int, N: Int, NRHS: Int,
                A: Array[Float], LDA: Int,
                B: Array[Float], LDB: Int,
                WORK: Array[Float], LWORK: Int, INFO: IntByReference)
      def dgels_(TRANS: String, M: Int, N: Int, NRHS: Int,
                A: Array[Double], LDA: Int,
                B: Array[Double], LDB: Int,
                WORK: Array[Double], LWORK: Int, INFO: IntByReference)
      def cgels_(TRANS: String, M: Int, N: Int, NRHS: Int,
                A: Array[Float], LDA: Int,
                B: Array[Float], LDB: Int,
                WORK: Array[Float], LWORK: Int, INFO: IntByReference)
      def zgels_(TRANS: String, M: Int, N: Int, NRHS: Int,
                A: Array[Double], LDA: Int,
                B: Array[Double], LDB: Int,
                WORK: Array[Double], LWORK: Int, INFO: IntByReference)

      def sgeev_(JOBVL: String, JOBVR: String,
                N: Int, A: Array[Float], LDA: Int,
                WR: Array[Float], WI: Array[Float],
                VL: Array[Float], LDVL: Int,
                VR: Array[Float], LDVR: Int,
                WORK: Array[Float], LWORK: Int, INFO: IntByReference)
      def dgeev_(JOBVL: String, JOBVR: String,
                N: Int, A: Array[Double], LDA: Int,
                WR: Array[Double], WI: Array[Double],
                VL: Array[Double], LDVL: Int,
                VR: Array[Double], LDVR: Int,
                WORK: Array[Double], LWORK: Int, INFO: IntByReference)
      def cgeev_(JOBVL: String, JOBVR: String,
                N: Int, A: Array[Float], LDA: Int,
                WR: Array[Float], WI: Array[Float],
                VL: Array[Float], LDVL: Int,
                VR: Array[Float], LDVR: Int,
                WORK: Array[Float], LWORK: Int, INFO: IntByReference)
      def zgeev_(JOBVL: String, JOBVR: String,
                N: Int, A: Array[Double], LDA: Int,
                WR: Array[Double], WI: Array[Double],
                VL: Array[Double], LDVL: Int,
                VR: Array[Double], LDVR: Int,
                WORK: Array[Double], LWORK: Int, INFO: IntByReference)
}
