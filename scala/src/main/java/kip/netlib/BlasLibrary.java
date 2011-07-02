package kip.netlib;

import com.sun.jna.Callback;
import com.sun.jna.Library;
import com.sun.jna.Native;
import com.sun.jna.NativeLibrary;
import com.sun.jna.NativeLong;
import com.sun.jna.Pointer;
import com.sun.jna.PointerType;
import com.sun.jna.Structure;
import com.sun.jna.ptr.DoubleByReference;
import com.sun.jna.ptr.FloatByReference;
import com.sun.jna.ptr.IntByReference;
import java.nio.DoubleBuffer;
import java.nio.FloatBuffer;
import java.nio.IntBuffer;


public interface BlasLibrary extends Library {
    public static final String JNA_LIBRARY_NAME = "vecLib";
    public static final BlasLibrary INSTANCE =
        (BlasLibrary)Native.loadLibrary(BlasLibrary.JNA_LIBRARY_NAME, BlasLibrary.class);
    
    // enum CBLAS_ORDER
    public static final int CblasRowMajor=101, CblasColMajor=102;
    
    // enum CBLAS_TRANSPOSE
    public static final int CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113, AtlasConj=114;
    
    public double cblas_dasum(int N, double[] X, int incX);

    public void cblas_dgemm(int Order, int TransA, int TransB,
                            int M, int N, int K,
                            double alpha,
                            double[] A, int lda,
                            double[] B, int ldb,
                            double beta,
                            double[] C, int ldc);

    public void cblas_cgemm(int Order, int TransA,
                            int TransB, int M, int N,
                            int K, double[] alpha, double[] A,
                            int lda, double[] B, int ldb,
                            double[] beta, double[] C, int ldc);
    
    //
    // LAPACK
    //
    
    
    public int dgels_(String trans,
                      IntByReference m, IntByReference n, 
                      IntByReference nrhs,
                      double[] a, IntByReference lda,
                      double[] b, IntByReference ldb, 
                      double[] work, IntByReference lwork, IntByReference info);

}
