package fftw3;

import com.sun.jna.*;
import java.nio.*;


public class Test {
    public static void main(String[] args) {
        
        FFTW3Library fft = FFTW3Library.INSTANCE;
        
        long sizeofComplex = 16;
        long N = 1024;
        long nbytes = sizeofComplex*N;
        Pointer in = fft.fftw_malloc(new NativeLong(nbytes));
        Pointer out = fft.fftw_malloc(new NativeLong(nbytes));
        DoubleBuffer inbuf = in.getByteBuffer(0, nbytes).asDoubleBuffer();
        DoubleBuffer outbuf = out.getByteBuffer(0, nbytes).asDoubleBuffer();
        
        int rank = 1;
        IntBuffer n = IntBuffer.wrap(new int[]{1024});
	FFTW3Library.fftw_plan p = fft.fftw_plan_dft(rank, n, inbuf, outbuf, FFTW3Library.FFTW_FORWARD, FFTW3Library.FFTW_ESTIMATE);
        fft.fftw_execute(p); 
        fft.fftw_destroy_plan(p);
        fft.fftw_free(in);
        fft.fftw_free(out);
    }
}
