package kip.md


import math._


trait Integrator {
  def dt: Double

/*  
  def run(time: Double): Double = {
    var i = 0
    while (i*dt < time) {
      singleStep()
      i += 1
    }
    i*dt
  }
  */
  
  def singleStep(world: World)
  
}


object Verlet {
  sealed abstract class Thermostat
  case object ThermoNone extends Thermostat
  case class ThermoLangevin(temp: Double, damp: Double, rand: util.Random) extends Thermostat
  case class ThermoRescaling(temp: Double) extends Thermostat
}

class Verlet(var dt: Double,
             var exclude: Set[Atom] = Set(),
             var thermostat: Verlet.Thermostat = Verlet.ThermoNone) extends Integrator {
  
  def singleStep(world: World) {
    for (a <- world.atoms) {
      a.x += dt*a.vx
      a.y += dt*a.vy
      a.z += dt*a.vz
    }
    world.volume.wrapAtoms(world.atoms)
    
    world.calculateForces()
    
    thermostat match {
      case Verlet.ThermoNone => ()
      case Verlet.ThermoLangevin(temp, damp, rand) => {
        for (a <- world.atoms) {
          val drag = a.mass / damp
          a.fx += - drag*a.vx + sqrt(2*temp*drag/dt)*rand.nextGaussian()
          a.fy += - drag*a.vy + sqrt(2*temp*drag/dt)*rand.nextGaussian()
          a.fz += - drag*a.vz + sqrt(2*temp*drag/dt)*rand.nextGaussian()
        }
      }
      case Verlet.ThermoRescaling(temp) => ()
    }
    
    for (a <- world.atoms) {
      val m = a.mass
      a.vx += dt*a.fx/m
      a.vy += dt*a.fy/m
      a.vz += dt*a.fz/m
    }
  }
}
