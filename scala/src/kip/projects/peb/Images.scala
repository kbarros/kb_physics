package kip.projects.peb

import kip.util.{Snapshot, LammpsParser}
import kip.util.Util.{formatDataInColumns, writeStringToFile}
import kip.util.Statistics.{mean}
import kip.math.Vec3
import kip.math.Math.sqr
import scala.math.sqrt



object Images {  
  def go(fileIn: String, fileOut: String, chainLength: Int, numChains: Int, time: Int, rCutoff: Double) {
    println("parsing "+fileIn+" at "+time)
    val snaps = LammpsParser.readLammpsDump(fileIn, maxSnapshots=1,
                                            process={s => if (s.time == time) Some(s) else None})
    val s = snaps(0)
    s.unwindCoords()
    val mols = PEB.molecules(chainLength, numChains)
    PEB.fixPolymerContinuity(s, mols)

    val monomerTypes = Set(1.0, 2.0) // graftedMonomer, monomer
    val counterionType = 3.0
    val saltCounterionType = 4.0
    val saltCoionType = 5.0
    
    val monomers        = s.typ.indices.filter { monomerTypes contains s.typ(_) }
    val counterions     = s.typ.indices.filter { counterionType == s.typ(_) }
    val saltCounterions = s.typ.indices.filter { saltCounterionType == s.typ(_) }
    val saltCoions      = s.typ.indices.filter { saltCoionType == s.typ(_) }
    
    def distanceFromTube(atom: Int): Double = {
      var minDist2 = java.lang.Double.MAX_VALUE
      for (mon <- monomers) {
        val dist2 = s.distance2(atom, mon)
        minDist2 = math.min(minDist2, dist2)
      }
      math.sqrt(minDist2)
    }

    var atomId = 1
    val sb = new StringBuilder()

    // print header
    val lx = s.hi.x-s.lo.x
    val ly = s.hi.y-s.lo.y
    val lz = s.hi.z-s.lo.z
    sb ++= "HEADER    POLYMERS\n"
    sb ++= "COMPND    LAMMPS\n"
    sb ++= "REMARK    Box: %5.1f x %5.1f x %5.1f\n".format(lx, ly, lz)
    
    def atomStr(idx: Int, label: String, p: Vec3): String = {
      val fmt = "ATOM  %5d  %s               %8.3f%8.3f%8.3f  1.00  0.00              \n"
      fmt.format(idx, label, p.x, p.y, p.z)
    }
    
    // print monomers
    for (i <- monomers) {
      sb ++= atomStr(atomId, "MN", s.principalImage(i))
      atomId += 1
    }

    // print monovalent counterions
    for (i <- counterions.filter(distanceFromTube(_) < rCutoff)) {
      sb ++= atomStr(atomId, "C1", s.principalImage(i))
      atomId += 1
    }

    // print salt counterions
    for (i <- saltCounterions.filter(distanceFromTube(_) < rCutoff)) {
      sb ++= atomStr(atomId, "CZ", s.principalImage(i))
      atomId += 1
    }

    // print salt coions
//    for (i <- saltCoions.filter(distanceFromTube(_) < rCutoff)) {
//      sb ++= atomStr(atomId, "CO", s.principalImage(i))
//      atomId += 1
//    }
    
    // print terminator
    sb ++= "END"
    
    kip.util.Util.writeStringToFile(sb.toString, fileOut)
  }
  
  def main(args: Array[String]) {
    go(fileIn="/home/kbarros/scratch/peb/run16.2x/z2Ns01Nm000r16-mino/dump.dat",
       fileOut="pdb/img1.pdb",
       chainLength=30, numChains=1, time=41000000, rCutoff=2)
    
    go(fileIn="/home/kbarros/scratch/peb/run16.2x/z2Ns15Nm000r16-mino/dump.dat",
       fileOut="pdb/img2.pdb",
       chainLength=30, numChains=1, time=12000000.toInt, rCutoff=2)
    
    go(fileIn="/home/kbarros/scratch/peb/run16.2x/z2Ns15Nm135r16-mino/dump.dat",
       fileOut="pdb/img3.pdb",
       chainLength=30, numChains=1, time=19000000.toInt.toInt, rCutoff=2)
  }
}
