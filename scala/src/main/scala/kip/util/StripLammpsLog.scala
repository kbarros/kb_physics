package kip.util

object StripLammpsLog {
  def main(args: Array[String]) {
    val (fname: String, minTime: Int) = args match {
      case Array(fname) => (fname, 0)
      case Array(fname, minTime) => (fname, minTime.toInt)
      case _ => println("Usage: StripLammpsLog filename [minTime]"); exit(1)
    }
    
    val thermo: Seq[Thermo] = LammpsParser.readLammpsThermo(fname)
    println("# time, temperature, potential, energy, pressure")
    thermo.foreach ( th =>
      if (th.time > minTime) {
        println(th.time+" "+th.temperature+" "+th.potential+" "+th.energy+" "+th.pressure)
      }
    )
  }
}
