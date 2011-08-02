package kip.md.apps


case class FrameData(index: Int, time: Double, temp: Double, natoms1: Int, force: Double, x: Array[Double], y: Array[Double], r: Array[Double], dislocations: (Array[Int], Array[Int]))

//object FrameData extends App {
//  val filename = "test.dat"
//  kip.util.Util.writeObjectGz(filename, FrameData(index=1, time=1.1, temp=1.2, natoms1=5, Array(1.0, 2.0), null, null, (null, null)))
//  val frame: FrameData = kip.util.Util.readObjectGz(filename)
//  println(frame)
//  println(frame.x.toList)
//}
