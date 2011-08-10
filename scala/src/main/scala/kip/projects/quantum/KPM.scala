package kip.projects.quantum

import smatrix._
import Constructors.complexDbl

// Kernel polynomial method



object KPM {
  def plotHistogram(a: Array[Double]) {
    val plot = new scikit.graphics.dim2.Plot("Histogram")
    scikit.util.Utilities.frame(plot.getComponent(), plot.getTitle())
    val hist = new scikit.dataset.Histogram(0.01)
    a.foreach { hist.accum(_) }
    plot.registerBars("Density", hist, java.awt.Color.BLACK)
  }
  
  // histogram eigenvalues  
  def eigenvaluesExact(m: PackedSparse[Scalar.ComplexDbl]): Array[Double] = {
    val (v, w) = m.toDense.eig
    v.map(_.re).toArray.sorted
  }
}

class KPM(m: PackedSparse[Scalar.ComplexDbl]) {
  
  // calculate moments
}
