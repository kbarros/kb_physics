package kip.math

import com.sun.jna.ptr.IntByReference


trait Blas {
  def nativeArray: Boolean
  
  def dgemm(TRANSA: String,TRANSB: String,
            M: Int, N: Int, K: Int,
            ALPHA: Double,
            A: Array[Double], LDA: Int,
            B: Array[Double], LDB: Int,
            BETA: Double,
            C: Array[Double], LDC: Int)
  
  def zgemm(TRANSA: String,TRANSB: String,
            M: Int, N: Int, K: Int,
            ALPHA: Array[Double],
            A: Array[Double], LDA: Int,
            B: Array[Double], LDB: Int,
            BETA: Array[Double],
            C: Array[Double], LDC: Int)
}

trait Lapack {
  def nativeArray: Boolean
  
  def dgels(TRANS: String, M: Int, N: Int, NRHS: Int,
            A: Array[Double], LDA: Int,
            B: Array[Double], LDB: Int,
            WORK: Array[Double], LWORK: Int, INFO: IntByReference)
  
  def zgels(TRANS: String, M: Int, N: Int, NRHS: Int,
            A: Array[Double], LDA: Int,
            B: Array[Double], LDB: Int,
            WORK: Array[Double], LWORK: Int, INFO: IntByReference)
}


object Netlib {
  import kip.nativelib.VecLib
  import com.sun.jna.Native
  
  object Implicits {
    implicit def intToIntByReference(a: Int) = new IntByReference(a)
  }
  
  lazy val blas: Blas = {
    val vecLib = Native.loadLibrary("vecLib", classOf[VecLib]).asInstanceOf[VecLib]
    
    new Blas {

      def transToEnum(t: String) = t match {
        case "N" | "n" => VecLib.CblasNoTrans
        case "T" | "t" => VecLib.CblasTrans
        case "C" | "c" => VecLib.CblasConjTrans
        case x => sys.error("Invalid matrix operation: "+x)
      }
      
      def nativeArray = true
      
      def dgemm(TRANSA: String, TRANSB: String,
            M: Int, N: Int, K: Int,
            ALPHA: Double,
            A: Array[Double], LDA: Int,
            B: Array[Double], LDB: Int,
            BETA: Double,
            C: Array[Double], LDC: Int) {
              vecLib.cblas_dgemm(VecLib.CblasColMajor, transToEnum(TRANSA), transToEnum(TRANSB),
                                 M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
      }
      
      def zgemm(TRANSA: String,TRANSB: String,
                M: Int, N: Int, K: Int,
                ALPHA: Array[Double],
                A: Array[Double], LDA: Int,
                B: Array[Double], LDB: Int,
                BETA: Array[Double],
                C: Array[Double], LDC: Int) {
        vecLib.cblas_zgemm(VecLib.CblasColMajor, transToEnum(TRANSA), transToEnum(TRANSB),
                           M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
      }
    }
  }

  lazy val lapack: Lapack = {
    val vecLib = Native.loadLibrary("vecLib", classOf[VecLib]).asInstanceOf[VecLib]
    new Lapack {
      
      import Implicits._
      
      def nativeArray = true
      
      def dgels(TRANS: String, M: Int, N: Int, NRHS: Int,
                A: Array[Double], LDA: Int,
                B: Array[Double], LDB: Int,
                WORK: Array[Double], LWORK: Int, INFO: IntByReference) {
                  vecLib.dgels_(TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO)
      }

      def zgels(TRANS: String, M: Int, N: Int, NRHS: Int,
                A: Array[Double], LDA: Int,
                B: Array[Double], LDB: Int,
                WORK: Array[Double], LWORK: Int, INFO: IntByReference) {
                  vecLib.zgels_(TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO)
      }
    }
  }

}
