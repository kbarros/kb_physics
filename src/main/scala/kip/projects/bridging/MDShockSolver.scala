package kip.projects.bridging


object MDStressFn extends StressFn {
  val md = new SubMD2d(ncols=2, nrows=2, a=1.12246, dt=0.01) // in production, we should use a 4x4 MD lattice or greater
  val rho0 = md.density0
  
  def stress(defgrad: Double, energy: Double): Double = {
    md.applyDeformationGradient(defgrad)
    md.initializeCrystalPositions()
    val zeroTemperature = true
    if (zeroTemperature) {
      md.applyKineticEnergy(0)
      md.convertCauchyToFirstPiolaKirchoff(md.virialStress())
    }
    else {
      md.applyKineticEnergy(energy*md.volume0 - md.potentialEnergy())
      val tLo = 1.0
      val tHi = 10.0
      var time = 0.0
      // skip transient behavior in (0 < t < tLo)
      while (time < tLo) {
        md.verletStep()
        time += md.dt
      }
      // average stress over times (tLo < t < tHi)
      // (testing is needed to set these parameters)
      var stressAcc = 0.0
      var stressAccCnt = 0
      while (time < tHi) {
        md.verletStep()
        stressAcc += md.convertCauchyToFirstPiolaKirchoff(md.virialStress())
        stressAccCnt += 1
        time += md.dt
      }
      stressAcc / stressAccCnt
    }
  }
  
  def zeroTempEnergy(defgrad: Double) = {
    md.applyDeformationGradient(defgrad)
    md.initializeCrystalPositions()
    md.applyKineticEnergy(0)
    md.potentialEnergy() / md.volume0
  }
}


class MDShockSolver(val L: Int, dx: Double, dt: Double, defgradAmplitude: Double, defgradWidth: Double) {
  // initial deformation gradient is a stretched region (L/4 < x < 3L/4) in a relaxed background
  import math.tanh
  val (a, w) = (defgradAmplitude, defgradWidth)
  def defgrad0 = Array.tabulate(L)(i => 1.0 + (a/2)*(tanh((dx/w)*(i-L/4)) - tanh((dx/w)*(i-3*L/4))))
  
  // two macro-solvers at first and second accuracy
  var s1 = new ElastodynamicLaws1d(L, MDStressFn.rho0, defgrad0, MDStressFn)  
  var s2 = new ElastodynamicLaws1d(L, MDStressFn.rho0, defgrad0, MDStressFn)  
  var time = 0.0
  
  def step() {
    CLIntegratorFirstOrder.iter(s1, dx, dt)
    CLIntegratorSecondOrder.iter(s2, dx, dt)
    time += dt
  }
}
