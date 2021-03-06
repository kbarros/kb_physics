package kip.projects.quantum.kpm2

import scala.reflect.ClassTag
import java.util.Arrays


class ArrayBuf[@specialized(Int, Float, Double) T: ClassTag]() {
  var size = 0
  var buffer = new Array[T](16)
  val scaleFactor = 1.5
  
  def resize(newSize: Int) {
    while (buffer.size < newSize) {
      val newSpace = (buffer.size*scaleFactor).round.toInt
      val newBuffer = new Array[T](newSpace)
      System.arraycopy(buffer, 0, newBuffer, 0, size)
      buffer = newBuffer
    }
    size = newSize
  }
  
  def push(x: T) {
    resize(size+1)
    buffer(size-1) = x
  }
  
  def clear() {
    size = 0
  }
  
  def apply(i: Int): T = {
    require(i < size)
    buffer(i)
  }
  
  def update(i: Int, x: T) {
    require(i < size)
    buffer(i) = x
  }
  
  def copyFrom(that: ArrayBuf[T]) {
    clear()
    resize(that.size)
    System.arraycopy(that.buffer, 0, buffer, 0, size)
  }
}


class DenseComplex(val numRows: Int, val numCols: Int, val data: Array[Double]) {
  val floatsPerElem = 2 // re + im
  require(data.size == floatsPerElem*numRows*numCols)
  def this(numRows: Int, numCols: Int) = this(numRows, numCols, new Array[Double](2*numRows*numCols))
  
  def zero() {
    for (i <- 0 until data.size)
      data(i) = 0.0
  }
  
  def dataIndex(i: Int, j: Int): Int = {
    j*numRows + i
  }
  
  def get_re(i: Int, j: Int): Double = {
    val k = dataIndex(i, j)
    data(floatsPerElem*k + 0)
  }

  def get_im(i: Int, j: Int): Double = {
    val k = dataIndex(i, j)
    data(floatsPerElem*k + 1)
  }

  def set(i: Int, j: Int, re: Double, im: Double) {
    val k = dataIndex(i, j)
    data(floatsPerElem*k + 0) = re
    data(floatsPerElem*k + 1) = im
  }
  
  def add(i: Int, j: Int, re: Double, im: Double) {
    val k = dataIndex(i, j)
    data(floatsPerElem*k + 0) += re
    data(floatsPerElem*k + 1) += im
  }
  
  // TODO: Remove!
  import smatrix.Dense
  import smatrix.Scalar
  import smatrix.Constructors
  import smatrix.Complexd
  def toSmatrix(): Dense[Scalar.ComplexDbl] = {
    val ret = Constructors.complexDbl.dense(numRows, numCols)
    for (i <- 0 until numRows;
         j <- 0 until numCols) {
      val k = dataIndex(i, j)
      ret(i, j) = Complexd(data(floatsPerElem*k+0), data(floatsPerElem*k+1))
    }
    ret
  }
  def fromSmatrix(that: Dense[Scalar.ComplexDbl]): this.type = {
    for (i <- 0 until numRows;
         j <- 0 until numCols) {
      val z = that(i, j)
      set(i, j, z.re, z.im)
    }
    this
  }
}


class SparseCooComplex(val numRows: Int, val numCols: Int) {
  var rowIdx = new ArrayBuf[Int]()
  var colIdx = new ArrayBuf[Int]()
  var data = new ArrayBuf[Double]()
    
  def clear() {
    rowIdx.clear()
    colIdx.clear()
    data.clear()
  }
  
  def add(i: Int, j: Int, re: Double, im: Double) {
    rowIdx.push(i)
    colIdx.push(j)
    data.push(re)
    data.push(im)
  }
  
  def toCsr(): SparseCsrComplex = {
    val ret = new SparseCsrComplex(numRows, numCols)
    ret.fromCoo(this)
    ret
  }
}


class SparseCsrComplex(val numRows: Int, val numCols: Int) {  
  // data sorted in row-major format
  var rowIdx = new ArrayBuf[Int]()
  var colIdx = new ArrayBuf[Int]()
  var data   = new ArrayBuf[Double]()
  val rowPtr = new Array[Int](numRows+1) 
  
  def nnz() = rowIdx.size
  
  def clear() {
    rowIdx.clear()
    colIdx.clear()
    data.clear()
  }
  
  def zero() {
    Arrays.fill(data.buffer, 0.0)
  }
  
  def fromCsr(that: SparseCsrComplex) {
    require(numRows == that.numRows && numCols == that.numCols)
    rowIdx.copyFrom(that.rowIdx)
    colIdx.copyFrom(that.colIdx)
    data.copyFrom(that.data)
    System.arraycopy(that.rowPtr, 0, rowPtr, 0, numRows+1)
  }
  
  def fromCoo(that: SparseCooComplex) {
    require(numRows == that.numRows && numCols == that.numCols)
    rowIdx.copyFrom(that.rowIdx)
    colIdx.copyFrom(that.colIdx)
    
    // sort indices
    def swap(a: Int, b: Int) {
      val ti = rowIdx(a)
      val tj = colIdx(a)
      rowIdx(a) = rowIdx(b)
      colIdx(a) = colIdx(b)
      rowIdx(b) = ti
      colIdx(b) = tj
    }
    def compare(a: Int, b: Int) = {
      val c1 = rowIdx(a) compare rowIdx(b)
      val c2 = colIdx(a) compare colIdx(b)
      if (c1 == 0) c2 else c1
    }
    Sort.quicksort(swap, compare, 0, nnz)
    
    // remove duplicate indices
    val numElemsWithDuplicates = that.rowIdx.size
    require(numElemsWithDuplicates > 0)
    var numUnique = 1
    for (i <- 1 until numElemsWithDuplicates) {
      if (compare(i-1, i) != 0) {
        rowIdx(numUnique) = rowIdx(i)
        colIdx(numUnique) = colIdx(i)
        numUnique += 1
      }
    }
    rowIdx.size = numUnique
    colIdx.size = numUnique
    
    // set row pointers
    var row = 0
    for (k <- 0 until numUnique) {
      while (row <= rowIdx(k)) {
        rowPtr(row) = k
        row += 1
      }
    }
    while (row <= numRows) {
      rowPtr(row) = numUnique
      row += 1
    }
    
    // fill data
    data.resize(2*numUnique)
    Arrays.fill(data.buffer, 0.0)
    for (k <- 0 until numElemsWithDuplicates) {
      val i = that.rowIdx(k)
      val j = that.colIdx(k)
      this += (i, j, that.data(2*k+0), that.data(2*k+1))
    }
  }
  
  def definedColumns(i: Int): IndexedSeq[Int] = {
    for (k <- rowPtr(i) until rowPtr(i+1)) yield colIdx(k)
  }
  
  def index(i: Int, j: Int): Int = {
    var k = rowPtr(i)
    // TODO: binary search
    for (k <- rowPtr(i) until rowPtr(i+1))
      if (j == colIdx(k))
        return k
    throw new java.lang.AssertionError(s"Undefined index ($i, $j)")
  }
  
  def set(i: Int, j: Int, re: Double, im: Double) {
    val k = index(i, j)
    data(2*k+0) = re
    data(2*k+1) = im
  }
  
  def get_re(i: Int, j: Int): Double = {
    val k = index(i, j)
    data(2*k+0)
  }
  
  def get_im(i: Int, j: Int): Double = {
    val k = index(i, j)
    data(2*k+1)
  }
  
  def get_abs2(i: Int, j: Int): Double = {
    val k = index(i, j)
    val re = data(2*k+0)
    val im = data(2*k+1)
    re*re + im*im
  }
  
  def trace_re(): Double = {
    require(numRows == numCols, "Trace operation requires square matrix")
    var acc = 0.0
    for (i <- 0 until numRows)
      acc += get_re(i, i)
    acc
  }
  
  def +=(i: Int, j: Int, re: Double, im: Double) {
    val k = index(i, j)
    data(2*k+0) += re
    data(2*k+1) += im
  }

  def +=(re: Double, im: Double) {
    require(numRows == numCols, "Add scalar operation requires square matrix")
    for (i <- 0 until numRows) {
      val k = index(i, i)
      data(2*k+0) += re
      data(2*k+1) += im
    }
  }
  
  def *=(a_re: Double, a_im: Double) {
    for (k <- 0 until nnz) {
      val d_re = data(2*k+0)
      val d_im = data(2*k+1)
      data(2*k+0) = d_re*a_re - d_im*a_im
      data(2*k+1) = d_re*a_im + d_im*a_re
    }
  }
  
  // TODO: Remove!
  import smatrix.PackedSparse
  import smatrix.Scalar
  import smatrix.Complexd
  def toSmatrix(): PackedSparse[Scalar.ComplexDbl] = {
    val ret = PackedSparse.fromIndices[Scalar.ComplexDbl](numRows, numCols, rowIdx.buffer.take(nnz) zip colIdx.buffer.take(nnz))
    for (k <- 0 until nnz) {
      ret(rowIdx(k), colIdx(k)) = Complexd(data(2*k+0), data(2*k+1))
    }
    ret
  }
  def fromSmatrix(that: PackedSparse[Scalar.ComplexDbl]): this.type = {
    val coo = new SparseCooComplex(that.numRows, that.numCols)
    for ((i, j) <- that.definedIndices) {
      val z = that(i, j)
      coo.add(i, j, z.re, z.im)
    }
    fromCoo(coo)
    this
  }
}


object MatrixOps {
  // C := alpha*A*B + beta*C
  def zcsrmm(alpha_re: Double, alpha_im: Double, A: SparseCsrComplex, B: DenseComplex, beta_re: Double, beta_im: Double, C: DenseComplex) {
    require(
        C.numRows == A.numRows &&
        A.numCols == B.numRows &&
        B.numCols == C.numCols, "Cannot multiply matrices of shape: [%d, %d] * [%d, %d] -> [%d, %d].".format(
            A.numRows, A.numCols, B.numRows, B.numCols, C.numRows, C.numCols))
    require(B ne C, "Illegal aliasing in matrix product.")
    for (i <- 0 until A.numRows) {
      if (beta_re != 1.0 || beta_im != 0.0) {
        for (j <- 0 until C.numCols) {
          val C_ij_re = C.get_re(i, j)
          val C_ij_im = C.get_im(i, j)
          val beta_C_ij_re = C_ij_re * beta_re - C_ij_im * beta_im
          val beta_C_ij_im = C_ij_re * beta_im + C_ij_im * beta_re
          C.set(i, j, beta_C_ij_re, beta_C_ij_im)
        }
      }
      
      var ptr = A.rowPtr(i)
      while (ptr < A.rowPtr(i+1)) {
        val k = A.colIdx(ptr)
        val A_ik_re = A.data(2*ptr + 0)
        val A_ik_im = A.data(2*ptr + 1)
        val alpha_A_ik_re = A_ik_re * alpha_re - A_ik_im * alpha_im
        val alpha_A_ik_im = A_ik_re * alpha_im + A_ik_im * alpha_re
        var j = 0
        while (j < C.numCols) {
          val idx = B.dataIndex(k, j)
          val B_kj_re = B.data(2*idx+0)
          val B_kj_im = B.data(2*idx+1)
          val re = alpha_A_ik_re * B_kj_re - alpha_A_ik_im * B_kj_im
          val im = alpha_A_ik_re * B_kj_im + alpha_A_ik_im * B_kj_re
          C.add(i, j, re, im)
          j += 1
        }
        ptr += 1
      }
    }
  }
}


object Sort {
  def quicksort(swap: (Int, Int) => Unit, compare: (Int, Int) => Int, off: Int, len: Int) {
    def sort1(l: Int, r: Int) {
      var pivot = (l + r) / 2
      var i = l; var j = r
      while (i <= j) {
        while (compare(i, pivot) < 0) i += 1
        while (compare(j, pivot) > 0) j -= 1
        if (i <= j) {
          swap(i, j)
          if      (pivot == i) pivot = j
          else if (pivot == j) pivot = i
          i += 1
          j -= 1
        }
      }
      if (l < j) sort1(l, j)
      if (j < r) sort1(i, r)
    }
    sort1(0, len - 1)
  }
}
