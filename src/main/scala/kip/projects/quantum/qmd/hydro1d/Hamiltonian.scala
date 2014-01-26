package kip.projects.quantum.qmd.hydro1d

import smatrix._
import kip.projects.quantum._
import kip.projects.quantum.ctor._
import math._
import kip.math.Vec3


// TODO:
//  - refactor; KPM gets own directory
//  - move random out of KPM
//  - implement Arpack, and move all energy scaling into KPM routines
//  - check energy derivatives
//  - try fourier transformation (figure out FFT mapping for BZ of non-square lattice)


object Hamiltonian {
  def main(args: Array[String]) {
    println("hello world")
    testEigenvalues()
  }
  
  def drawMatrix(m: Dense[S]) {
    val n = m.numRows
    val data = Array.tabulate[Double](n*n) { k =>
      val i = (n-1) - k / n
      val j = k % n
      m(i, j).abs
    }
    scikit.util.Commands.grid(n, n, data)
  }
  
  def mapMatrix(m: Dense[S], f: R => R): Dense[S] = {
    val (w, v) = m.eig
    val mp = v * diag(w.map(x => f(x.re): Complexf)) * v.tran
    mp
  }
  
  def testEigenvalues() {
    val n = 50
    
    val h = (new Hamiltonian(numAtoms=n, decay=2, cutoff=4.0, rand=new util.Random(0))).matrix()
    
    val mu = 0.3
    val fn_energy:  (R => R) = e => if (e < mu) (e - mu) else 0
    val fn_filling: (R => R) = e => if (e < mu) 1.0 else 0
    
    drawMatrix(h.toDense)
    drawMatrix(mapMatrix(h.toDense, fn_filling))
    drawMatrix(mapMatrix(h.toDense, fn_energy))
    
    val order = 1000
    val range = KPM.range(npts=5*order)
    
    val eig = KPM.eigenvaluesExact(h)
    val exactDensity = KPM.integrateDeltas(range, eig, moment=0).map(_ / n)
    
    val kpm = new KPM(h, nrand=10, seed=0)
    val r = kpm.randomVector()
    val moments = kpm.momentsStochastic(order, r)
    val kpmEigenvalues = range.map(KPM.densityOfStates(moments, KPM.jacksonKernel(order), _))
    val kpmDensity = KPM.integrate(range, kpmEigenvalues, moment=0).map(_ / n)
    
    val plot = KPM.mkPlot("Integrated density of states")
    KPM.plotLines(plot, (range, exactDensity), "Exact", java.awt.Color.RED)
    KPM.plotLines(plot, (range, kpmDensity), "Approx", java.awt.Color.BLACK)
    
    val exactEnergy  = eig.map(fn_energy(_)).sum
    val exactFilling = eig.map(fn_filling(_)).sum
    val c_energy  = KPM.expansionCoefficients2(order, quadPts=10*order, fn_energy)
    val c_filling = KPM.expansionCoefficients2(order, quadPts=10*order, fn_filling)
    val kpmEnergy  = (c_energy,  moments).zipped.map(_*_).sum
    val kpmFilling = (c_filling, moments).zipped.map(_*_).sum
    println(s"kpmEnergy  = $kpmEnergy ($exactEnergy)")
    println(s"kpmFilling = $kpmFilling ($exactFilling)")
  }
}

class Hamiltonian(val numAtoms: Int, val decay: Double, val cutoff: Double, val rand: util.Random) {
  val pos = new Array[Double](numAtoms)
  
  randomizePositions(spacing=1.0, disorder=0.0)
  
  def randomizePositions(spacing: Double, disorder: Double) {
    for (i <- 0 until numAtoms) yield {
      pos(i) = spacing*i + disorder*rand.nextGaussian()
    }
  }
  
  def grouping(i: Int, s: Int): Int = {
    require(numAtoms / s == 0)
    i % s
  }
  
  def matrix(): PackedSparse[S] = {
    val ret = sparse(numAtoms, numAtoms): HashSparse[S]
    for (i <- 0 until numAtoms;
         j <- 0 until numAtoms) {
      val r = abs(pos(i) - pos(j))
      if (r < cutoff) {
        val t_ij = exp(-r/decay)
        ret(i, j) = t_ij
        ret(j, i) = t_ij
      }
    }
    
    // TODO: use Arpack && move this into KPM routines
    val eig = KPM.eigenvaluesExact(ret.toPacked)
    val e_max = eig.max
    val e_min = eig.min
    val e_avg   = (e_max + e_min)/2
    val e_scale = (e_max - e_min)/2
    for (i <- 0 until numAtoms) { ret(i,i) -= e_avg }
    ret /= (1.01*e_scale)
    
    ret.toPacked
  }

  // returns gradient (dE/dr_i)
  def fieldDerivative(dEdH: PackedSparse[S]): Array[R] = {
    null
  }
}
