package kip.util

import java.io._

object DCDParser {  
  case class Frame(lx: Double, ly: Double, lz: Double, x: Array[Double], y: Array[Double], z: Array[Double])
  
  val bytes4 = new Array[Byte](4)
  val bytes8 = new Array[Byte](8)
  
  def readInt(in: InputStream): Int = {
    val b = bytes4
    in.read(b)
    ((b(3)&0xff)<<24) + ((b(2)&0xff)<<16) + ((b(1)&0xff)<<8) + (b(0)&0xff)
  }
  
  def readLong(in: InputStream): Long = {
    val i1 = readInt(in)
    val i2 = readInt(in)
    (i2.toLong << 32) + i1.toLong
  }
  
  def readFloat(in: InputStream): Float = {
    java.lang.Float.intBitsToFloat(readInt(in))
  }

  def readDouble(in: InputStream): Double = {
    java.lang.Double.longBitsToDouble(readLong(in))
  }
  
  def parseFile(file: String): Array[Frame] = {
    val in = new BufferedInputStream(new FileInputStream(file))
    
    // Read file header
    assert(readInt(in) == 84) // first magic no.
    assert(readInt(in) == 1146244931) // second magic no. ("CORD" in ascii)
    val numFrames = readInt(in)
    val startTimestep = readInt(in) // initial timestep
    val skipTimesteps = readInt(in) // timesteps between frames in file
    val numTimesteps = readInt(in) // number of timesteps in simulation
    assert(readInt(in) == 0)
    assert(readInt(in) == 0)
    assert(readInt(in) == 0)
    assert(readInt(in) == 0)
    assert(readInt(in) == 0)
    assert(readInt(in) == 0) // timestep (unused)
    assert(readInt(in) == 1) // include unit cell
    assert(readInt(in) == 0)
    assert(readInt(in) == 0)
    assert(readInt(in) == 0)
    assert(readInt(in) == 0)
    assert(readInt(in) == 0)
    assert(readInt(in) == 0)
    assert(readInt(in) == 0)
    assert(readInt(in) == 0)
    assert(readInt(in) == 24) // Pretend to be CHARMM version 24
    assert(readInt(in) == 84) // ?
    assert(readInt(in) == 164) // ?
    assert(readInt(in) == 2) // ?
    val b = new Array[Byte](80)
    in.read(b)
    println(b.map(_.toChar).mkString)
    in.read(b)
    println(b.map(_.toChar).mkString)
    assert(readInt(in) == 164)
    assert(readInt(in) == 4)
    val numParticles = readInt(in)
    assert(readInt(in) == 4)
    
    
    def readFrame(): Frame = {
      // frame header
      assert(readInt(in) == 48)
      val lx = readDouble(in)
      assert(readDouble(in) == 0f) // 90 degrees
      val ly = readDouble(in)
      assert(readDouble(in) == 0f) // 90 degrees
      assert(readDouble(in) == 0f) // 90 degrees
      val lz = readDouble(in)
      assert(readInt(in) == 48)
      
      // frame data
      val Seq(x, y, z) = for (i <- 0 until 3) yield {
        val FLOAT_BYTES = 4
        assert(readInt(in) == numParticles*FLOAT_BYTES)
        val x = Array.fill(numParticles) { readFloat(in).toDouble }
        assert(readInt(in) == numParticles*FLOAT_BYTES)
        x
      }
      Frame(lx, ly, lz, x, y, z)
    }
    
    val frames = Util.time("Reading data")(Array.fill(numFrames)(readFrame()))
    in.close()
    frames
  }
  
  def visualizer(frames: Seq[Frame]) = {
    import kip.md._
    import kip.graphics._
    import kip.math._
    
    val f = frames(0)
    
    val viz = new Visualizer()
    viz.setBounds(Bounds3d(Vec3(-f.lx/2,-f.ly/2,-f.lz/2), Vec3(f.lx/2, f.ly/2, f.lz/2)))

    new {
      def goto(i: Int) {
        val f = frames(i)
        viz.setString("Frame "+i)
        val spheres = for (j <- 0 until f.x.size) yield {
          // assert (f.x(j) >= 0)
          Visualizer.Sphere(Vec3(f.x(j), f.y(j), f.z(j)), radius=0.5, color=java.awt.Color.BLUE)
        }
        viz.setParticles(spheres)
        viz.display()
      }
    }
  }
}
