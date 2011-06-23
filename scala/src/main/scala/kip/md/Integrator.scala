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
      a.pos.x += dt*a.v.x
      a.pos.y += dt*a.v.y
      a.pos.z += dt*a.v.z
    }
    world.volume.wrapAtoms(world.atoms)
    
    world.calculateForces()
    
    thermostat match {
      case Verlet.ThermoNone => ()
      case Verlet.ThermoLangevin(temp, damp, rand) => {
        for (a <- world.atoms) {
          val drag = a.mass / damp
          a.f.x += - drag*a.v.x + sqrt(2*temp*drag/dt)*rand.nextGaussian()
          a.f.y += - drag*a.v.y + sqrt(2*temp*drag/dt)*rand.nextGaussian()
          a.f.z += - drag*a.v.z + sqrt(2*temp*drag/dt)*rand.nextGaussian()
        }
      }
      case Verlet.ThermoRescaling(temp) => ()
    }
    
    for (a <- world.atoms) {
      val m = a.mass
      a.v.x += dt*a.f.x/m
      a.v.y += dt*a.f.y/m
      a.v.z += dt*a.f.z/m
    }
  }
}
