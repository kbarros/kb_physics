package kip.projects.quantum.kpm2

import scala.reflect.ClassTag


class ArrayBuf[@specialized(Int, Float, Double) T: ClassTag]() {
  var size = 0
  var buffer = new Array[T](16)
  val scaleFactor = 1.5
  
  def grow(len: Int) {
    while (buffer.size < len) {
      val newSpace = (buffer.size*scaleFactor).round.toInt
      val newBuffer = new Array[T](newSpace)
      System.arraycopy(buffer, 0, newBuffer, 0, size)
      buffer = newBuffer
    }
  }
  
  def add(x: T) {
    grow(size+1)
    buffer(size) = x
    size += 1
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
    grow(that.size)
    System.arraycopy(that.buffer, 0, buffer, 0, that.size)
    size = that.size
  }
}


class DenseComplex(val numRows: Int, val numCols: Int, val data: Array[Double]) {
  val floatsPerElem = 2 // re + im
  require(data.size == floatsPerElem*numRows*numCols)
  def this(numRows: Int, numCols: Int) = this(numRows, numCols, new Array[Double](2*numRows*numCols))
  
  def clear() {
    for (i <- 0 until data.size)
      data(i) = 0.0
  }
  
  def dataIndex(i: Int, j: Int): Int = {
    j*numRows + i
  }
  
  def set(i: Int, j: Int, re: Double, im: Double) {
    val k = dataIndex(i, j)
    data(floatsPerElem*k + 0) = re
    data(floatsPerElem*k + 1) = im
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
  
  def nnz() = rowIdx.size
  
  def clear() {
    rowIdx.clear()
    colIdx.clear()
    data.clear()
  }
  
  def add(i: Int, j: Int, re: Double, im: Double) {
    rowIdx.add(i)
    colIdx.add(j)
    data.add(re)
    data.add(im)
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
    data.copyFrom(that.data)
    
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
    
    // set row pointers
    var row = 0
    for (k <- 0 until nnz) {
      while (row <= rowIdx(k)) {
        rowPtr(row) = k
        row += 1
      }
    }
    while (row <= numRows) {
      rowPtr(row) = nnz
      row += 1
    }
    
    // fill data
    for (k <- 0 until nnz) {
      val i = that.rowIdx(k)
      val j = that.colIdx(k)
      set(i, j, that.data(2*k+0), that.data(2*k+1))
    }
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
  
  def +=(a_re: Double, a_im: Double) {
    require(numRows == numCols, "Add scalar operation requires square matrix")
    for (i <- 0 until numRows) {
      val k = index(i, i)
      data(2*k+0) += a_re
      data(2*k+1) += a_im
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
          if (pivot == i) pivot = j
          if (pivot == j) pivot = i
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
