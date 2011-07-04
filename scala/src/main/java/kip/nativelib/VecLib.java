package kip.nativelib;

import com.sun.jna.Library;
import com.sun.jna.ptr.IntByReference;


public interface VecLib extends Library {
    // enum CBLAS_ORDER
    public static final int CblasRowMajor=101, CblasColMajor=102;
    
    // enum CBLAS_TRANSPOSE
    public static final int CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113, AtlasConj=114;
    
    //
    // BLAS
    //
    
    public double cblas_dasum(int N, double[] X, int incX);

    public void cblas_dgemm(int Order, int TransA, int TransB,
                            int M, int N, int K,
                            double alpha,
                            double[] A, int lda,
                            double[] B, int ldb,
                            double beta,
                            double[] C, int ldc);

    public void cblas_zgemm(int Order, int TransA, int TransB,
                            int M, int N, int K,
                            double[] alpha,
                            double[] A, int lda,
                            double[] B, int ldb,
                            double[] beta,
                            double[] C, int ldc);
    
    //
    // LAPACK
    //
    
    public int dgels_(String trans,
                      IntByReference m, IntByReference n, 
                      IntByReference nrhs,
                      double[] a, IntByReference lda,
                      double[] b, IntByReference ldb, 
                      double[] work, IntByReference lwork, IntByReference info);

    public int zgels_(String trans,
                      IntByReference m, IntByReference n, 
                      IntByReference nrhs,
                      double[] a, IntByReference lda,
                      double[] b, IntByReference ldb, 
                      double[] work, IntByReference lwork, IntByReference info);

    public int dgeev_(String jobvl, String jobvr, IntByReference n,
                      double[] a, IntByReference lda,
                      double[] wr, double[] wi,
                      double[] vl, IntByReference ldvl,
                      double[] vr, IntByReference ldvr,
                      double[] work, IntByReference lwork, IntByReference info);

    public int zgeev_(String jobvl, String jobvr, IntByReference n,
                      double[] a, IntByReference lda,
                      double[] wr, double[] wi,
                      double[] vl, IntByReference ldvl,
                      double[] vr, IntByReference ldvr,
                      double[] work, IntByReference lwork, IntByReference info);
}
