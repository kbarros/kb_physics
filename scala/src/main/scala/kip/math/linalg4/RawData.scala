package kip.math.linalg4
import java.nio.FloatBuffer
import java.nio.DoubleBuffer
import java.nio.Buffer


object RawData {
  
  object Builder {
    implicit val FltArray = new Builder[Float, FloatBuffer] {
      def build(size: Int) = new FltArray(size)
    }
    implicit val DblArray = new Builder[Double, DoubleBuffer] {
      def build(size: Int) = new DblArray(size)
    }
  }
  
  trait Builder[@specialized(Float, Double) Raw, Buf <: Buffer] {
    def build(size: Int): RawData[Raw, Buf]
  }  


  class FltArray(size: Int) extends RawData[Float, FloatBuffer] {
    val array = new Array[Float](size)
    val buffer = java.nio.FloatBuffer.wrap(array)
    def apply(i: Int): Float = array(i)
    def update(i: Int, x: Float) { array(i) = x }
  }

  class DblArray(size: Int) extends RawData[Double, DoubleBuffer] {
    val array = new Array[Double](size)
    val buffer = java.nio.DoubleBuffer.wrap(array)
    def apply(i: Int): Double = array(i)
    def update(i: Int, x: Double): Unit = array(i) = x
  }
}


trait RawData[@specialized(Float, Double) Raw, Buf <: Buffer] {
  def buffer: Buf
  def apply(i: Int): Raw
  def update(i: Int, x: Raw)
  def dispose() = ()
}

