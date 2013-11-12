package kip.projects.magnet_dynamics

import kip.graphics.Bounds3d
import kip.graphics.RetainedScene
import kip.math.Vec3
import scikit.graphics.dim2.Grid
import scikit.jobs.Control
import scikit.jobs.Job
import scikit.jobs.Simulation
import kip.projects.kondo_rkky.SpinViz

object Spirals extends App {
  new Control(new Spirals(), "Skirmions")  
}

class Spirals extends Simulation {
  var lattice: TriLattice = _
  var sim: MagnetDynamics = _
  val rs = new RetainedScene(Bounds3d(Vec3(0,0,0),Vec3(1,1,0)), sizew=600, sizeh=600, cameraDistance=1.2)
      
  // we work with magnetic moments |S|=1, but physical system
  // has spin 3/2 hbar. this becomes important when converting
  // external magnetic field from Tesla to meV units.
  val electronMagneticMoment = 5.7883818066e-2 // meV/T
  val siteMagneticMoment = 3 * electronMagneticMoment  // meV/T
  
  
  def load(c: Control) {
    params.add("Lx", 21)
    params.add("Ly", 21)
    params.add("Lz", 1)
    
    // couplings in units meV (with magnetic moment |S|=1) 
    params.addm("J1 (meV)", -2.8)
    params.addm("J2 (meV)", -0.48)
    params.addm("J3 (meV)", -0.08)
    params.addm("Jz (meV)", 0.02)
    params.addm("Dx (meV)", -0.59)
    params.addm("Dz (meV)", 0.48)
    params.addm("Hz (T)", 0)
    
//    params.addm("J1 (meV)", -2.8)
//    params.addm("J2 (meV)", 0.0)
//    params.addm("J3 (meV)", 0.0)
//    params.addm("Jz (meV)", 0.0)
//    params.addm("Dx (meV)", 0.0)
//    params.addm("Dz (meV)", 0.0)
//    params.addm("Hz (T)", 0.0)
    
    params.addm("T (meV)", 0.01)
    params.addm("alpha", 0.1)
    params.addm("dt", 1e-1)
    
    params.add("energy")
    params.add("net spin")
  }

  def animate() {
    lattice.J1 = params.fget("J1 (meV)")
    lattice.J2 = params.fget("J2 (meV)")
    lattice.J3 = params.fget("J3 (meV)")
    lattice.Jz = params.fget("Jz (meV)")
    lattice.Dx = params.fget("Dx (meV)")
    lattice.Dz = params.fget("Dz (meV)")
    lattice.Hz = params.fget("Hz (T)") * siteMagneticMoment
    
    sim.T = params.fget("T (meV)")
    sim.alpha = params.fget("alpha")
    sim.dt = params.fget("dt")
    
    val field = new Array[Double](3*lattice.N)
    for (i <- 0 until lattice.N) {
      field(0 + 3*i) = sim.s.x(i)
      field(1 + 3*i) = sim.s.y(i)
      field(2 + 3*i) = sim.s.z(i)
    }
    val viz = new SpinViz(lattice.Lx, lattice.Ly)
    viz.drawSpins(field, rs)
    
    var acc_x = sim.s.x.sum / sim.s.x.length
    var acc_y = sim.s.y.sum / sim.s.y.length
    var acc_z = sim.s.z.sum / sim.s.z.length
    
    params.set("net spin", math.sqrt(acc_x*acc_x + acc_y*acc_y + acc_z*acc_z))
    
    params.set("energy", lattice.energy(sim.s))
  }

  def clear() {
  }

  def run() {
    val Lx = params.iget("Lx")
    val Ly = params.iget("Ly")
    val Lz = params.iget("Lz")

    lattice = new TriLattice(
        Lx=Lx, Ly=Ly, Lz=Lz,
        J1=params.fget("J1 (meV)"),
        J2=params.fget("J2 (meV)"),
        J3=params.fget("J3 (meV)"),
        Jz=params.fget("Jz (meV)"),
        Dx=params.fget("Dx (meV)"),
        Dz=params.fget("Dz (meV)"),
        Hz=params.fget("Hz (T)") * siteMagneticMoment)
    
    sim = new MagnetDynamics(
        T=params.fget("T (meV)"),
        alpha=params.fget("alpha"),
        dt=params.fget("dt"),
        lattice=lattice)
    
    while (true) {
      Job.animate()
      
      for (i <- 0 until 5)
        // sim.stepEuler()
        sim.stepHeun()
    }
  }
}
