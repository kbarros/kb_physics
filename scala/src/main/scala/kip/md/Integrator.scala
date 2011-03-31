package kip.md


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


class Verlet(var dt: Double, var exclude: Set[Atom] = Set()) extends Integrator {
  
  def singleStep(world: World) {
    for (a <- world.atoms) {
      a.x += dt*a.vx
      a.y += dt*a.vy
      a.z += dt*a.vz
    }
    world.volume.wrapAtoms(world.atoms)
    
    world.calculateForces()
    
    for (a <- world.atoms) {
      val m = a.mass
      a.vx += dt*a.fx/m
      a.vy += dt*a.fy/m
      a.vz += dt*a.fz/m
    }
  }
}
