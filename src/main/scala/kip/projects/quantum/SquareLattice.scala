package kip.projects.quantum

import smatrix._
import ctor._
import kip.math.Vec3


object SquareLattice extends App {
  import kip.util.Util.{time}
  // time("integrated density")(testIntegratedDensity())
  time("spiral energy")(testSpiralEnergy())
  
  // Plots the integrated density of states
  def testIntegratedDensity() {
    val q = new SquareLattice(w=20, h=20, t1=1, t2= -0.5, J_H=0.1, e_min= -10, e_max= 10)
    q.setFieldFerro(q.field)
    q.fillMatrix(q.matrix)
    
    val H = q.matrix
    require((H - H.dag).norm2.abs < 1e-10, "Found non-hermitian hamiltonian!")
    println("N = "+H.numRows)
    val order = 1000
    val kpm = new KPM(H, nrand=1, seed=0)
    val range = KPM.range(npts=5*order)
    
    val plot = KPM.mkPlot("Integrated density of states")
    KPM.plotLines(plot, (range, KPM.integrateDeltas(range, KPM.eigenvaluesExact(H), moment=0)), "Exact", java.awt.Color.RED)
    KPM.plotLines(plot, (range, KPM.integrate(range, KPM.eigenvaluesApprox(order, range, kpm), moment=0)), "Approx", java.awt.Color.BLACK)
  }
  
  def testSpiralEnergy() {
    val mu = -2.5
    
    def energyPerSite(q: KondoHamiltonian) = {
      val eigs = KPM.eigenvaluesExact(q.matrix).map(_ * q.e_scale) 
      eigs.takeWhile(_ <= mu).map(_.re - mu).sum / q.numLatticeSites
    }

    def testJ(J: Double) {
      println(s"Using J=$J")
      
      println("d=6")
      for (L <- Seq(12, 18, 24, 30, 36)) {
        val q = new SquareLattice(w=L, h=L, t1=1, t2= -0.5, J_H=J, e_min= -10, e_max=10)
        q.setFieldSpiral(q.field, 6)
        q.fillMatrix(q.matrix)
        println(s"L=$L, e=${energyPerSite(q)}")
      }

      println("d=8")
      for (L <- Seq(8, 16, 24, 32, 40)) {
        val q = new SquareLattice(w=L, h=L, t1=1, t2= -0.5, J_H=J, e_min= -10, e_max=10)
        q.setFieldSpiral(q.field, 8)
        q.fillMatrix(q.matrix)
        println(s"L=$L, e=${energyPerSite(q)}")
      }
    }
    
    // testJ(0.1)
    // testJ(0.05)
    // testJ(0)
  }
}

class SquareLattice(val w: Int, val h: Int, val t1: R, val t2: R, val J_H: R, val e_min: R, val e_max: R) extends KondoHamiltonian {
  require(h % 2 == 0, "Need even number of rows, or hamiltonian is non-hermitian")
  
  val numLatticeSites = w*h
  
  def coord2idx(x: Int, y: Int): Int = {
    x + y*w
  }
  
  def idx2coord(i: Int) = {
    (i%w, i/w)
  }
  
  def setField(desc: String) {
    desc match {
      case "ferro"    => setFieldFerro(field)
      case s if desc.startsWith("spiral") => {
        setFieldSpiral(field, desc.replace("spiral", "").toInt)
      }
    }
  }
    
  def setFieldSpiral(field: Array[R], d: Int) {
    require(w % d == 0)
    require(h % d == 0)
    for (y <- 0 until h;
         x <- 0 until w) {
      val i = coord2idx(x, y)      
      field(fieldIndex(0, i)) = math.cos(2*math.Pi * (x+y) / d)
      field(fieldIndex(1, i)) = math.sin(2*math.Pi * (x+y) / d)
      field(fieldIndex(2, i)) = 0
    }
  }

  
  lazy val latticePositions: Array[Vec3] = {
    (for (i <- 0 until numLatticeSites) yield { 
      val (x, y) = idx2coord(i)
      Vec3(x, y, 0)
    }).toArray
  }

  // returns (x, y) indices for the `nn`th neighbor site
  def neighbor1(x: Int, y: Int, nn: Int): (Int, Int) = {
    val delta = Array((1, 0), (0, 1), (-1, 0), (0, -1))
    val (dx, dy) = delta(nn)
    ((x+dx+w)%w, (y+dy+h)%h)
  }

  // returns (x, y) indices for the `nn`th neighbor site
  def neighbor2(x: Int, y: Int, nn: Int): (Int, Int) = {
    val delta = Array((2, 0), (0, 2), (-2, 0), (0, -2))
    val (dx, dy) = delta(nn)
    ((x+dx+w)%w, (y+dy+h)%h)
  }
  
  // in *unscaled* (physical) energy units 
  val hoppingMatrix: PackedSparse[S] = {
    val ret = sparse(2*numLatticeSites, 2*numLatticeSites)
    ret.clear()
    
    for (y <- 0 until h;
         x <- 0 until w;
         nn <- 0 until 4;
         sp <- 0 until 2) {
      val i = matrixIndex(sp, coord2idx(x, y))
      
      val (nx1, ny1) = neighbor1(x, y, nn);
      val j1 = matrixIndex(sp, coord2idx(nx1, ny1))
      ret(i, j1) = -t1
      
      val (nx2, ny2) = neighbor2(x, y, nn);
      val j2 = matrixIndex(sp, coord2idx(nx2, ny2))
      ret(i, j2) = -t2
    }
    ret.toPacked
  }

}
