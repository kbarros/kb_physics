package kip.util

import java.io._

object QHull {
  
  case class Vertices2d(a: Array[Double]) {
    val n: Int = a.size/2
    def x(i: Int): Double = a(2*i+0)
    def y(i: Int): Double = a(2*i+1)
  }
  
  def delaunay2d(vertices: Vertices2d): Array[Int] = {
    val p = Runtime.getRuntime().exec("qdelaunay QJ i Fx")
    val pw = new PrintWriter(p.getOutputStream())
    val br = new BufferedReader(new InputStreamReader(p.getInputStream()));
    
    // write vertices to qhull
    pw.println(2) // dimension
    pw.println(vertices.n)
    for (i <- 0 until vertices.n) {
      pw.println(vertices.x(i)+" "+vertices.y(i))
    }
    pw.close()
    
    // read triangle output from "i" option
    val ntris = br.readLine().toInt
    val tris = new Array[Int](3*ntris)
    for (i <- 0 until ntris) {
      val Array(v1, v2, v3) = br.readLine().split("\\s+")
      tris(3*i+0) = v1.toInt
      tris(3*i+1) = v2.toInt
      tris(3*i+2) = v3.toInt
    }
    
    tris
  }
  

  def test() = {
    val a = Array[Double](0,0, 0,1, 1,1, 1,0, 0.6, 0.7)
    delaunay2d(Vertices2d(a))
  }
}
