package kip.projects.network

import scala.collection.mutable
import scala.util.Random
import scala.math._

object Network {
  
  type Rank = Int // MPI rank
  type Node = Int // index of node in cluster
  type Trace = Map[Rank, Seq[Msg]] // messages associated with each rank
  type Assignment = Array[Node] // assignment of ranks to nodes
  
  case class Msg(r1: Rank, r2: Rank, size: Double)

  case class Topology(w: Int, h: Int, d: Int, periodic: Boolean = true) {
    val size = w*h*d // number of nodes
    
    def node2coords(node: Node): (Int, Int, Int) = {
      val n = node
      val x = n % w
      val y = (n / w) % h
      val z = (n / w / h) % d
      (x, y, z)
    }
    
    def coords2node(x: Int, y: Int, z:Int): Node = {
      z*w*h + y*w + x
    }
    
    def delta(node1: Node, node2: Node): (Int, Int, Int) = {
      val c1 = node2coords(node1)
      val c2 = node2coords(node2)
      val dx = c2._1 - c1._1
      val dy = c2._2 - c1._2
      val dz = c2._3 - c1._3
      def wrap(d: Int, len: Int) = if (d > len/2) (d - len) else if (d < -len/2) (d + len) else d
      if (periodic)
        (wrap(dx, w), wrap(dy, h), wrap(dz, d))
      else
        (dx, dy, dz)
    }
    
    def hops(node1: Node, node2: Node): Double = {
      val (dx, dy, dz) = delta(node1, node2)
      abs(dx) + abs(dy) + abs(dz)
    }
    
    def format(fn: Node => String, colWidth: Int = 8): String = {
      val bldr = new StringBuilder
      def writeColumn(s: String) {
        bldr ++= s
        bldr ++= (Seq.fill(colWidth-s.size)(' '))
      }
      for (z <- 0 until d) {
        for (y <- 0 until h) {
          for (x <- 0 until w) {
            val node = coords2node(x, y, z)
            writeColumn(fn(node))
          }
          bldr += '\n'
        }
        bldr += '\n'
      }
      bldr.toString
    }
  }
  
  case class Config(top: Topology, trace: Trace, assign: Assignment) {
    def cost(rank: Rank): Double = {
      val msgs = trace(rank)
      (for (m <- msgs) yield {
        val n1 = assign(m.r1)
        val n2 = assign(m.r2)
        m.size * top.hops(n1, n2)
      }).sum
    }
    
    def cost: Double = {
      assign.indices.map(cost(_)).sum
    }
    
    def formatAssignment: String = {
      val node2rank = assign.zipWithIndex.toMap
      top.format(n => node2rank(n).toString)
    }
    
    def formatCost: String = {
      val node2rank = assign.zipWithIndex.toMap
      top.format(n => cost(node2rank(n)).toString)
    }
  }
  
  
  def swapElems[T](i: Int, j: Int, array: Array[T]) {
    val temp = array(i)
    array(i) = array(j)
    array(j) = temp
  }
    

  def readTrace(file: String): Trace = {
    val trace = new mutable.HashMap[Rank, mutable.ArrayBuffer[Msg]]()
    for (line <- scala.io.Source.fromFile(file).getLines()) {
      val tokens = line.trim.split("""\s+""")
      val msg: Option[Msg] = tokens.size match {
        case 0 => None
        case 2 => Some(Msg(tokens(0).toInt, tokens(1).toInt, 1d))
        case 3 => Some(Msg(tokens(0).toInt, tokens(1).toInt, tokens(2).toDouble))
        case _ => assert(false); None
      }
      def insertMsg(r: Rank, m: Msg) {
        if (!trace.contains(r))
          trace(r) = mutable.ArrayBuffer()
        trace(r) += m
      }
      msg.foreach { m =>
        insertMsg(m.r1, m)
        insertMsg(m.r2, m)
      }
    }
    
    trace.toMap.mapValues(_.toSeq)
  }
  
  def randomAssignment(numRanks: Int, rand: Random): Assignment = {
    val assign: Assignment = Array.tabulate(numRanks)(identity)
    for (i <- assign.indices) {
      val j = i + rand.nextInt(assign.size-i)
      swapElems(i, j, assign)
    }
    assign
  }
}


object MCOptimize {
  import Network._
  
  val rand = new Random(1)
  val top = Topology(w=4, h=4, d=4)
  val trace = readTrace(sys.props("user.home")+"/Desktop/network/4x4x4-mesh.dat")
  val assign = randomAssignment(numRanks=trace.keys.max+1, rand)
  val conf = Config(top, trace, assign)
  
  val plot = new scikit.graphics.dim2.Plot("Convergence")
  val frame = scikit.util.Utilities.frame(plot.getComponent(), plot.getTitle());
  val history = new scikit.dataset.DynamicArray()
  
  def metropolisStep(beta: Double, conf: Config, rand: Random) {
    val r1 = rand.nextInt(conf.assign.size)
    val r2 = rand.nextInt(conf.assign.size)
    
    val energy1 = conf.cost(r1) + conf.cost(r2)
    swapElems(r1, r2, conf.assign)
    val energy2 = conf.cost(r1) + conf.cost(r2)
    
    val de = energy2 - energy1
    if (de < 0 || exp(-beta*de) > rand.nextDouble()) {
      // leave rank assignments swapped
    }
    else {
      // revert swap
      swapElems(r1, r2, conf.assign)
    }
  }
  
  def go() {
    val betas = Seq(0.5, 0.8, 0.9, 1.0, 1.1, 1.2)
    var t = 0.0
    for (beta <- betas) {
      for (i <- 0 until 500) {
        for (j <- 0 until 1000) {
          metropolisStep(beta, conf, rand)
        }
        t += 1
        history.append2(t, conf.cost)
        plot.registerPoints("Data", history, java.awt.Color.BLACK)
      }
    }
  }
}
