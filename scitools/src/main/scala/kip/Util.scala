package kip

import java.lang.Math._
import java.io.{BufferedWriter, FileWriter}


object Util {
    final def sqr(x: Double) =
	x*x
    
    final def printTime[A](fn: => A) = {
	val begin = System.currentTimeMillis
	print ("Begin timing... ")
	val ret = fn
	val end = System.currentTimeMillis
	println("done. Execution time: " + (end-begin)/1000. + "s")
	ret
    }
    
    def average(vs: Seq[Double]) =
	vs.reduceLeft(_+_) / vs.length
    
    def trace[A](str: String, v: A) = {
	println(str + v)
	v
    }
    
    def fft1d_continuous(a: Array[Double], L:Double) = {
	import scikit.numerics.fft.FFT1D
	import scikit.dataset.Accumulator

	val accumulator = new Accumulator()
	
	val fft = new FFT1D(a.length)
	fft.setLength(L)
	fft.transform(a, new FFT1D.MapFn() {
	    def apply(k:Double, re:Double, im:Double) {
		accumulator.accum(k, re*re+im*im)
	    }
	})
	accumulator
    }
    
    def formatDataInColumns(kvs: (String, Array[Double])*) = {
	val sb = new StringBuffer()

	val (descs, vals) = List.unzip(kvs)
	sb.append("# " + descs.mkString(" ") + "\n")
	for (i <- 0 until vals(0).size) {
	    sb.append(vals.map{_(i)}.mkString(" ")+"\n")
	}
	sb.toString()
    }

    def writeStringToFile(s: String, fn: String) {
	val writer = new BufferedWriter(new FileWriter(fn))
    	writer.write(s);
    	writer.close();
    }
    
    def plot(a: Array[Double]) {
	scikit.util.Commands.plot(a)
    }
}
