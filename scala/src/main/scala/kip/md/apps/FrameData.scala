package kip.md.apps

import net.liftweb.json
import java.io.FileOutputStream
import java.util.zip.ZipOutputStream
import java.util.zip.ZipEntry
import java.io.BufferedWriter
import java.io.OutputStreamWriter
import java.util.zip.GZIPOutputStream
import java.io.FileInputStream
import java.util.zip.GZIPInputStream
import java.io.BufferedReader
import java.io.InputStreamReader


object FrameData {
  case class Frame(index: Int, time: Double, x: Array[Double], y: Array[Double], r: Array[Double], dislocations: (Array[Int], Array[Int]))
  
  class Writer(filename: String) {
    implicit val formats = json.Serialization.formats(json.NoTypeHints)
    val fout = new GZIPOutputStream(new FileOutputStream(filename))
    val writer = new BufferedWriter(new OutputStreamWriter(fout));

    def writeFrame(f: Frame) {
      json.Serialization.write(f, writer)
    }

    def close() {
      writer.close();
    }
  }
 
  class Reader(filename: String) {
    implicit val formats = json.Serialization.formats(json.NoTypeHints)
    val fin = new GZIPInputStream(new FileInputStream(filename))
    val reader = new BufferedReader(new InputStreamReader(fin));

    def readFrame(): Frame = {
      json.Serialization.read[Frame](reader)
    }

    def close {
      reader.close();
    }
  }

}


object FrameDataTest extends App {
  val reader = new FrameData.Reader(args(0))
  
  while (true) {
    println(reader.readFrame().index)
  }
}