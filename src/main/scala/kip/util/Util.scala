package kip.util

import scala.util.matching.Regex
import java.io.{LineNumberReader, BufferedReader, BufferedWriter, FileWriter, FileReader}
import java.io.FileOutputStream
import java.io.ObjectOutputStream
import java.io.FileInputStream
import java.io.ObjectInputStream
import java.util.zip.GZIPOutputStream
import java.util.zip.GZIPInputStream
import java.io.File
import java.io.IOException


object Util {
  def time[A](s: String)(f: => A): A = {
    print("Timing '"+s+"'...")
    val t1 = System.currentTimeMillis
    val ret = f
    val t2 = System.currentTimeMillis
    println(" done. Elapsed time "+(t2-t1)/1000.+"s");
    ret
  }

  def notime[A](s: String)(f: => A): A = {
    f
  }

  def tr[A](v: A, str: String) = {
    println(str + v)
    v
  }

  def formatDataInColumns(kvs: (String, Traversable[Double])*) = {
    val sb = new StringBuffer()

    val (descs, cols) = kvs.unzip
    val rows = cols.transpose
    sb.append("# " + descs.mkString(" ") + "\n")
    for (row <- rows) {
      sb.append(row.mkString(" ")+"\n")
    }
    sb.toString()
  }

  def readStringFromFile(fn: String) = {
    val stream = new java.io.FileInputStream(new java.io.File(fn))
    try {
      val fc = stream.getChannel()
      val bb = fc.map(java.nio.channels.FileChannel.MapMode.READ_ONLY, 0, fc.size())
      java.nio.charset.Charset.defaultCharset().decode(bb).toString()
    }
    finally {
      stream.close();
    }
  }
  
  def writeStringToFile(s: String, fn: String) {
    val writer = new BufferedWriter(new FileWriter(fn))
    writer.write(s);
    writer.close();
  }
  
  def readCsvFile(fn: String): (Array[Array[Double]], String) = {
    import scala.collection.mutable.ArrayBuffer
    var data: Array[ArrayBuffer[Double]] = null
    var desc: String = ""
    
    val reader = new LineNumberReader(new FileReader(fn))
    var line = reader.readLine()
    while (line != null) {
      line = line.trim
      
      if (line.size == 0) {
        // skip
      }
      else if (line.startsWith("#")) {
        desc += line
      }
      else {
        val items = line.split("\\s+").map(_.toDouble)
        if (data == null) {
          data = Array.fill(items.size)(new ArrayBuffer[Double]())
        }
        if (data.size != items.size) {
          throw new Error("Column size mismatch, "+data.size+"!="+items.size+", in '"+fn+"':"+reader.getLineNumber())
        }
        for (i <- items.indices) {
          data(i).append(items(i))
        }
      }
      line = reader.readLine()
    }
    
    (data.map(_.toArray), desc) 
  }

  def writeObjectGz[A <: Serializable](filename: String, f: A) {
    val out = new ObjectOutputStream(new GZIPOutputStream(new FileOutputStream(filename)))
    out.writeObject(f);
    out.close();
  }

  def readObjectGz[A <: Serializable](filename: String): A = {
    val in = new ObjectInputStream(new GZIPInputStream(new FileInputStream(filename)))
    val ret = in.readObject().asInstanceOf[A]
    in.close()
    ret
  }

  def createEmptyDir(dir: File) {
    if (dir.exists) {
      if (!dir.isDirectory()) {
        throw new IOException("File '%s' exists but is not directory".format(dir))
      }
      if (dir.list.size > 0) {
        println("WARNING: Using non-empty directory " + dir)
      }
    }
    if (!dir.exists) {
      dir.mkdir()
    }
  }
}
