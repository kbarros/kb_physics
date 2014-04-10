package kip.projects.quantum.qmd

import kip.util.Util.time
import kip.util.Interpreter
import kip.math.{Vec3}
import kip.graphics._
import kip.enrich._
import java.awt.{BorderLayout, Color, Frame, FileDialog}
import java.awt.event.{WindowAdapter, WindowEvent}
import javax.swing.{JPanel, JSlider, JComponent}
import javax.swing.event.{ChangeEvent, ChangeListener}
import java.io.File
import scikit.graphics.dim2.Plot
import kip.projects.quantum.kpm.KPMUtil
import kip.projects.quantum.kpm.EnergyScale
import scikit.dataset.PointSet
import kip.util.RangeArray
import kip.util.Snapshot


object TbViz {
  def readSnap(file: File): TbMD.Snap = {
    kip.util.JacksonWrapper.deserialize[TbMD.Snap](file.slurp)
  }
  
  def readSnaps(dir: File, readEvery: Int): Array[TbMD.Snap] = {
    require(dir.isDirectory(), "Cannot load directory %s".format(dir))  
    val files = dir.listFiles()
    val indices = files.map { f =>
      f.getName().takeWhile(_ != '.').toInt
    }
    (files zip indices).sortBy(_._2).filter(_._2 % readEvery == 0).map(x => readSnap(x._1))
  }
  
  def convertToLammpsSnapshot(snap: TbMD.Snap): Snapshot = {
    val natoms = snap.x.size / 3
    val ret = new Snapshot(time = 0, natoms)
    ret.lo = Vec3(snap.bdsLo(0), snap.bdsLo(1), snap.bdsLo(2))
    ret.hi = Vec3(snap.bdsHi(0), snap.bdsHi(1), snap.bdsHi(2))
    ret.id  = Array.tabulate(natoms)(i=>i.toDouble)
    ret.typ = Array.fill(natoms)(0.0)
    val pos = snap.x.grouped(3).toArray
    ret.x = pos.map(_(0)).toArray
    ret.y = pos.map(_(1)).toArray
    ret.z = pos.map(_(2)).toArray
    // ix, iy, iz
    val vel = snap.v.grouped(3)
    ret.vx = vel.map(_(0)).toArray
    ret.vy = vel.map(_(1)).toArray
    ret.vz = vel.map(_(2)).toArray
    // q
    ret
  }
  
  def main(args: Array[String]) {
    if (args.size != 2) {
      println("Usage: TbViz <dirname> <readevery>")
    }
    val dir = args(0)
    val readEvery = args(1).toInt
    val snaps = readSnaps(new java.io.File(dir+"/dump"), readEvery)
    val snaps2 = snaps.map(convertToLammpsSnapshot(_))
    val mv = new MolViz(snaps2, RenderProperties.basic)
    // new TbViz(snaps)
    
    val dv = new DensityViz(snaps)
    val pv = new PairViz(snaps2)
    mv.callbacks :+= (dv.display(_)) 
    mv.callbacks :+= (pv.display(_))
    mv.goto(0)
  }
}

class DensityViz(snaps: Array[TbMD.Snap]) {
  val plot = {
    val ret = new Plot("Integrated density")
    scikit.util.Utilities.frame(ret.getComponent(), ret.getTitle())
    ret
  }
  
  def display(idx: Int) {
    val snap = snaps(idx)
    val es = new EnergyScale(snap.energyScale._1, snap.energyScale._2)
    val gamma = KPMUtil.momentTransform(snap.moments, 10*snap.moments.size)
    val (xp, irho) = KPMUtil.integratedDensityFunction(gamma, es)
    val data = new PointSet(xp, irho)
    plot.registerLines("Integrated density", data, Color.BLACK)
    val muIdx = xp.indexWhere(_ > snap.mu)
    val muData = new PointSet(Array(snap.mu), Array(irho(muIdx)))
    plot.registerBars("Mu", muData, Color.RED)
  }
}


class PairViz(snaps: Array[Snapshot]) {
  def pairCorrelation(snaps: Seq[Snapshot], dr: Double, rmax: Double, typs1: Seq[Int], typs2: Seq[Int]): RangeArray[Double] = {
    val volume = snaps(0).volume
    val g = RangeArray.fill(xmin=0, xmax=rmax, dx=dr)(0d)
    // sum over all snapshots, and all pairs of particles. distances are binned into g
    var uniquePairs = 0
    var atomsCnt = 0
    
    val sameTyps = typs1.toSet == typs2.toSet
    
    for ((s,iter) <- snaps.zipWithIndex) {
      val ids1 = (0 until s.natoms) filter (i => typs1.contains(s.typ(i)))
      val ids2 = (0 until s.natoms) filter (i => typs2.contains(s.typ(i)))
      // if (iter % 100 == 0) println("Processing snapshot "+iter)
      for (i1 <- ids1; i2 <- ids2; if i1 != i2) {
        val r = math.sqrt(s.distance2(i1, i2))
        if (g.isDefinedAt(r))
          g(r) += 1
        uniquePairs += 1
      }
      atomsCnt += (ids1 ++ ids2).distinct.size
    }
    if (false) {
      for (r <- g.elemCenters) {
        val volume_fraction = 4*math.Pi*r*r*dr / volume
        g(r) /= volume_fraction * uniquePairs
      }
    }
    else {
      var acc = 0.0
      for (r <- g.elemCenters) {
        acc += g(r) / atomsCnt
        g(r) = acc
      }
    }
    g
  }
  
  
  val plot = {
    val ret = new Plot("Pair correlation")
    scikit.util.Utilities.frame(ret.getComponent(), ret.getTitle())
    ret
  }

  def display(idx: Int) {
    val snap = snaps(idx)
    val g = pairCorrelation(Seq(snap), dr=0.005, rmax=8, typs1=Seq(0), typs2=Seq(0))
    val data = new PointSet(g.elemCenters, g.elems)
    plot.registerLines("Integrated density", data, Color.BLACK)
  }
}
