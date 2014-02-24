package kip.projects.quantum.qmd

object Units {
  // length
  val angstrom = 1.0
  val bohr = 0.5292*angstrom // 0.5291772109217

  // energy
  val eV = 1.0
  val rydberg = 13.605804*eV // 13.6056925330
  
  // mass
  val amu = 1.0
  val massSi = 28.0851*amu
  
  // time
  val fs = math.sqrt(amu/eV) * angstrom / 10.1805054836529
  
  // temperature
  val kelvin = 1.0
  val kB = 8.6173323849609e-5 * eV / kelvin
}
