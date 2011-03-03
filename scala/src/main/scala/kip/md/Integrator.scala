package kip.md


trait Integrator {
  var dt: Double
  def step()
  def time: Double
}


class Verlet(world: World, var dt: Double) extends Integrator {
  var time = 0.0
  
  def step() {
    for (a <- world.atoms) {
      a.x += dt*a.vx
      a.y += dt*a.vy
      a.z += dt*a.vz
    }
    world.wrapAtoms()
    
    world.calculateForces()
    
    for (a <- world.atoms) {
      val m = a.mass
      a.vx += dt*a.fx/m
      a.vy += dt*a.fy/m
      a.vz += dt*a.fz/m
    }
    
    time += dt
  }
}
