package kip.projects.quantum.qmd.hydro1d

import smatrix._
import kip.projects.quantum._
import kip.projects.quantum.ctor._
import math._
import kip.math.Vec3

// TODO:
//  - check energy derivatives

object Hamiltonian {
  def main(args: Array[String]) {
    testEigenvalues()
  }
  
  def drawMatrix(mat: Dense[S]) {
    val n = mat.numRows
    val m = mat.numCols
    val data = Array.tabulate[Double](n*m) { k =>
      val i = (n-1) - k / m
      val j = k % m
      mat(i, j).abs
    }
    scikit.util.Commands.grid(m, n, data)
  }
  
  def mapMatrix(m: Dense[S], f: R => R): Dense[S] = {
    val (w, v) = m.eig
    val mp = v * diag(w.map(x => f(x.re): Complexf)) * v.tran
    mp
  }

  def fourierMatrix(n: Int): Dense[S] = {
    val ret = dense(n, n)
    ret.tabulate { (i, j) =>
      ((I*2*Pi*i*j)/n).exp / sqrt(n)
    } 
    ret
  }
  
  def testEigenvalues() {
    val lx = 36
    val ly = 36
    val n = lx*ly
    val order = 1000
    val decay = 2.0
    val cutoff = 4.0
    val spacing = 1.0
    val disorder = 0.2
    val mu = 0.3
    val nrand = 64
    
    val hamiltonian = new Hamiltonian(lx=lx, ly=ly, decay=decay, cutoff=cutoff, rand=new util.Random(0))
    hamiltonian.randomizePositions(spacing=spacing, disorder=disorder)
    val h = hamiltonian.matrix()
    
    val fn_energy:  (R => R) = e => if (e < mu) (e - mu) else 0
    val fn_filling: (R => R) = e => if (e < mu) 1.0 else 0
    
    val energyMatrix = mapMatrix(h.toDense, fn_energy)
    val densityMatrix = mapMatrix(h.toDense, fn_filling)
    val exactEnergy  = energyMatrix.trace
    val exactFilling = densityMatrix.trace
    
    drawMatrix(h.toDense)
    drawMatrix(densityMatrix)
    drawMatrix(energyMatrix)
    
    val kpm = new KPM(h, nrand=nrand, seed=2)
    val c_energy  = KPM.expansionCoefficients2(order, quadPts=10*order, fn_energy)
    val c_filling = KPM.expansionCoefficients2(order, quadPts=10*order, fn_filling)
    
    def estimateErrorsEmpirical(genRand: () => Dense[S]) {
      val iters = 1000
      val (es, ns) = (for (i <- 0 until iters) yield {
        val r = genRand()
        val moments = kpm.momentsStochastic(order, r)
        val energy  = (c_energy,  moments).zipped.map(_*_).sum
        val filling = (c_filling, moments).zipped.map(_*_).sum
        (energy.toDouble, filling.toDouble)
      }).unzip
      import kip.util.Statistics._
      println(s"energy  ${mean(es)}+-${stddev(es)} ($exactEnergy)")
      println(s"filling ${mean(ns)}+-${stddev(ns)} ($exactFilling)")
    }
//    println("correlated vectors")
//    estimateErrorsEmpirical(() => kpm.correlatedVectors(hamiltonian.grouping))
//    println("uncorrelated")
//    estimateErrorsEmpirical(() => kpm.randomVector())
    
    def estimateErrorsDeterministic(nr: Int) {
      val r  = dense(n, nr)
      val a = Array[S#A](1, I, -1, -I).map(_ * math.sqrt(nr))
      r.tabulate { (i, j) =>
        if (hamiltonian.grouping(i, nr) == j) a(kpm.rand.nextInt(4)) else 0
      }
      // val r = kpm.correlatedVectors(hamiltonian.grouping)
      
      val q = (r * r.dag) / nr
      
      var acc_e_c = 0.0
      var acc_e_uc = 0.0
      var acc_n_c = 0.0
      var acc_n_uc = 0.0
      for ((i, j) <- densityMatrix.definedIndices) {
        if (i != j) {
          acc_e_c  += (energyMatrix(i,j) * q(i, j)).abs2
          acc_e_uc += (energyMatrix(i,j) / sqrt(nr)).abs2
          acc_n_c  += (densityMatrix(i,j) * q(i, j)).abs2
          acc_n_uc += (densityMatrix(i,j) / sqrt(nr)).abs2
        }
      }
      // println(s"e = $exactEnergy +- (c ${sqrt(acc_e_c)}) (uc ${sqrt(acc_e_uc)})")
      // println(s"n = $exactFilling +- (c ${sqrt(acc_n_c)}) (uc ${sqrt(acc_n_uc)})")
      println(s"$nr ${sqrt(acc_e_c)} ${sqrt(acc_e_uc)} ${sqrt(acc_n_c)} ${sqrt(acc_n_uc)}")
    }
    println("# nrand sig_e_c sig_e_uc sig_n_c sig_n_uc")
    for (nr <- Seq(1, 4, 9, 16, 36, 81, 144, 324)) {
      estimateErrorsDeterministic(nr)
    }
    
    val r = kpm.correlatedVectors(hamiltonian.grouping)
    drawMatrix(r * r.dag)
    
//    val eig = KPM.eigenvaluesExact(h)
//    val moments = kpm.momentsStochastic(order, r)
//    val range = KPM.range(npts=5*order)
//    val kpmEigenvalues = range.map(KPM.densityOfStates(moments, KPM.jacksonKernel(order), _))
//    val kpmDensity = KPM.integrate(range, kpmEigenvalues, moment=0).map(_ / n)        
//    val exactDensity = KPM.integrateDeltas(range, eig, moment=0).map(_ / n)
//    val plot = KPM.mkPlot("Integrated density of states")
//    KPM.plotLines(plot, (range, exactDensity), "Exact", java.awt.Color.RED)
//    KPM.plotLines(plot, (range, kpmDensity), "Approx", java.awt.Color.BLACK)
  }
}

class Hamiltonian(val lx: Int, val ly: Int, val decay: Double, val cutoff: Double, val rand: util.Random) {
  val numAtoms = lx*ly
  val pos = new Array[Vec3](numAtoms)
  
  def latticeCoords(i: Int) = (i%lx, i/lx)
  
  def randomizePositions(spacing: Double, disorder: Double) {
    for (i <- 0 until numAtoms) {
      val (x, y) = latticeCoords(i)
      val dx = disorder*rand.nextGaussian()
      val dy = if (ly == 1) 0 else disorder*rand.nextGaussian()
      pos(i) = Vec3(x+0.5*y+dx, y*sqrt(3.0)/2.0+dy, 0)*spacing
    }
  }
  
  def grouping(i: Int, s: Int): Int = {
    if (ly == 1) {
      require(numAtoms % s == 0)
      i % s
    }
    else {
      val len = sqrt(s).toInt
      require(len*len == s)
      val (x, y) = latticeCoords(i)
      val xp = x % len
      val yp = y % len
      xp + yp*len
    }
  }
  
  def matrix(): PackedSparse[S] = {
    val ret = sparse(numAtoms, numAtoms): HashSparse[S]
    for (i <- 0 until numAtoms;
         j <- 0 until numAtoms) {
      val delta = pos(i) - pos(j)
      val dx = min(abs(delta.x), abs(lx-delta.x))
      val dy = min(abs(delta.y), abs(ly*sqrt(3.0)/2.0-delta.y))
      val r = sqrt(dx*dx + dy*dy)
      // val r = (pos(i) - pos(j)).norm
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
