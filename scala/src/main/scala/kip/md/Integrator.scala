package kip.md


trait Integrator {
  def step(dt: Double)
}


class Verlet(world: World) extends Integrator {
  
  def step(dt: Double) {
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
