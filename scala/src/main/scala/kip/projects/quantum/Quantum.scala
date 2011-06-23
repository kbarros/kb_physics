package kip.projects.quantum

import scalala.scalar._
import scalala.tensor.::
import scalala.tensor.mutable._
import scalala.tensor.dense._
import scalala.tensor.sparse._
import scalala.operators.Implicits._

import scalala.collection.sparse._

import scala.collection.mutable._


class SparseVectorC(rank: Int) {
  val elemsRe = new SparseArray[Double](rank)
  val elemsIm = new SparseArray[Double](rank)
  
  def *(that: DenseVectorC): Complex = {
    
    val activeCols = elemsRe.indexArray
    val reArray = elemsRe.valueArray
    val imArray = elemsIm.valueArray
    
    var re_acc = 0d
    var im_acc = 0d
    for ((col,i) <- activeCols.zipWithIndex) {
      val re1 = reArray(i)
      val im1 = imArray(i)
      val re2 = that.elemsRe(col)
      val im2 = that.elemsIm(col)
      
      re_acc += re1*re2 - im1*im2
      im_acc += re1*im2 + re2*im1
    }
    
    Complex(re_acc, im_acc)
  }
  
  def update(i: Int, c: Complex) {
    elemsRe(i) = c.real
    elemsIm(i) = c.imag
  }
  
  def apply(i: Int): Complex = {
    Complex(elemsRe(i), elemsIm(i))
  }
}

class SparseMatrixC(rank: Int) {
  val rows = Array.fill(rank)(new SparseVectorC(rank))
  
  def update(i: Int, j: Int, c: Complex): Unit = {
    rows(i)(j) = c
  }
  
  def apply(i: Int, j: Int): Complex = {
    rows(i)(j)
  }
  
  def *(v: DenseVectorC): DenseVectorC = {
    val ret = new DenseVectorC(rank)
    for (i <- 0 until rank) {
      ret(i) = rows(i) * v
    }
    ret
  }
  
  def *(that: DenseMatrixC): DenseMatrixC = {
    val ret = new DenseMatrixC(rank)
    for (j <- 0 until rank) {
      ret(::, j) := (this * that(::, j))
    }
    ret
  }
}


class DenseVectorC(val rank: Int, val elemsRe: DenseVectorCol[Double], val elemsIm: DenseVectorCol[Double]) {
  def this(rank: Int) {
    this(rank, DenseVectorCol.zeros[Double](rank), DenseVectorCol.zeros[Double](rank))
  }
  
  def update(i: Int, c: Complex) {
    elemsRe(i) = c.real
    elemsIm(i) = c.imag
  }
  
  def apply(i: Int): Complex = {
    Complex(elemsRe(i), elemsIm(i))
  }
  
  def :=(that: DenseVectorC) {
    elemsRe := that.elemsRe
    elemsIm := that.elemsIm
  }
}

class DenseMatrixC(val rank: Int, val elemsRe: DenseMatrix[Double], val elemsIm: DenseMatrix[Double]) {
  def this(rank: Int) {
    this(rank, DenseMatrix.zeros[Double](rank, rank), DenseMatrix.zeros[Double](rank, rank))
  }
  
  def apply(i: scalala.tensor.SelectAll, j: Int): DenseVectorC = {
    new DenseVectorC(rank, elemsRe(::,  j), elemsIm(::, j))
  }
}


trait CalculateTrace {
  def calc(f: DenseVectorC => DenseVectorC): Double
}

object StochasticTrace extends CalculateTrace {
  def calc(f: DenseVectorC => DenseVectorC): Double = {
    0
  }
}


// file format:
// (row col val)
// 

object Quantum {
  
}
