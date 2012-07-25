package kip.projects.bridging

object StandaloneMDShockApp extends App {
  var sim: MDShockSolver = _
  val L = 400
  val dx = 1.0
  val dt = 0.02
  sim = new MDShockSolver(L=L, dx=dx, dt=dt, defgradAmplitude=0.04, defgradWidth=5*dx)
  for (i <- 0 until 50) {
    sim.step()
    
    // general outputs include momentum, energy, and deformation gradient fields.
    // here, we simply print out the energy at position L/4
    // (according to first and second order macro-scale solvers)
     println("time: %g, energy1 %g, energy2 %g".format(sim.time, sim.s1.energy(L/4), sim.s2.energy(L/4)))
  }
}
