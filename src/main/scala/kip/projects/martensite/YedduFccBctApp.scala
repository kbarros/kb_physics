package kip.projects.martensite

import scikit.graphics.dim2.Grid
import scikit.jobs.Control
import scikit.jobs.Job
import scikit.jobs.Simulation
import scikit.jobs.params.DoubleValue
import scikit.jobs.params.Parameters
import scikit.numerics.Math2
import scala.util.Random
import scala.math.sqrt
import scikit.graphics.dim2.Plot
import scikit.dataset.PointSet


// 3d vector
case class vec3(var x: Double, var y: Double, var z: Double) {
    def *(c: Double) = vec3(c*x, c*y, c*z)
}

// symmetric tensor
case class symten(var xx: Double, var yy: Double, var zz: Double, var xy: Double, var yz: Double, var zx: Double) {
  def *(c: Double) = symten(c*xx, c*yy, c*zz, c*xy, c*yz, c*zx)
}


class FccBctSim(params: Parameters) {
  var timeSteps = 0
  
  // Convert everything to SI units (meter, second, kg, joule, pascal)
  val L = params.fget("L (um)") * 1e-6
  val L_nuc = params.fget("L_nuc (um)") * 1e-6
  val n0 = params.iget("n0")
  val nx = n0
  val ny = n0
  val nz = n0
  val dx = L / n0
  val Dg = params.fget("Dg (J/m^3)")                  // "driving force"
  val gamma = params.fget("gamma (J/m^2)")            // interfacial energy
  val delta = params.fget("delta / dx") * dx          // interfacial width
  val E = params.fget("E (GPa)") * 1e9                // youngs modulus
  val nu = params.fget("nu")                          // poisson ratio
  val b3 = params.fget("b3")                          // bain strains
  val b1 = params.fget("b1")
  val rho = params.fget("rho (kg/m^3)")               // material density
  val v = math.sqrt(E*(1-nu) / (rho*(1+nu)*(1-2*nu))) // compression wave speed
  val tau = L / v                                     // elastic propagation time
  val dt = params.fget("dt * (v / dx)") * dx / v
  val sigma_y = params.fget("sigma_y (MPa)") * 1e6    // yield stress
  val gamma0 = params.fget("gamma0 / tau") * tau      // damping time
  val M = params.fget("L delta Dg / v") * v / (delta * math.abs(Dg)) // mobility
  val k = params.fget("k / tau") * tau
  
  val N = nx*ny*nz
  
  val W = 3*gamma/delta
  println(s"W=$W dg=$Dg")
  val beta = 6*gamma*delta
  val A_hat = W
  val B_hat = 2*W - 4*Dg
  val C_hat = W - 3*Dg
  
  val random = new Random(params.iget("Random seed", 0))

  // displacement vector
  val u = Array.fill(N)(vec3(0, 0, 0))
  
  // du/dt
  val udot = Array.fill(N)(vec3(0, 0, 0))
  
  // khachaturyan fields eta1, eta2, eta3
  val eta = Array(Array.fill(N)(0.0),
                  Array.fill(N)(0.0),
                  Array.fill(N)(0.0))
  
  // bain strains
  val bain = Array(symten(b3, b1, b1, 0, 0, 0),
                   symten(b1, b3, b1, 0, 0, 0),
                   symten(b1, b1, b3, 0, 0, 0))
  
  // plastic strain
  val eps_pl = Array.fill(N)(symten(0, 0, 0, 0, 0, 0))
  
  // these fields derive from the above dynamical ones
  val eps_tot = Array.fill(N)(symten(0, 0, 0, 0, 0, 0)) // total strain
  val eps_tr  = Array.fill(N)(symten(0, 0, 0, 0, 0, 0)) // transformation strain
  val eps_el  = Array.fill(N)(symten(0, 0, 0, 0, 0, 0)) // elastic strain
  val sigma   = Array.fill(N)(symten(0, 0, 0, 0, 0, 0)) // stress
  val divsig  = Array.fill(N)(vec3(0, 0, 0))            // divergence of stress
  val laplace_eta = Array(Array.fill(N)(0.0),           // laplacian of eta
                          Array.fill(N)(0.0),
                          Array.fill(N)(0.0))
  val el_energy = Array.fill(N)(0.0)
  val gl_energy = Array.fill(N)(0.0)
  
  
  // utility lattice functions
  def wrap(x: Int, n: Int) = ((x % n) + n) % n
  def idxToCoords(i: Int) = ((i % nx), ((i / nx) % ny), (i / (nx*ny)))
  def coordsToIdx(x: Int, y: Int, z: Int) = wrap(x, nx) + wrap(y, ny)*nx + wrap(z, nz)*nx*ny
  
  // arrays containing nearest neighbor indices
  def mkOffsetArray(dx: Int, dy: Int, dz: Int) = Array.tabulate(N) { i: Int =>
    val (x,y,z) = idxToCoords(i)
    assert(x < nx)
    assert(y < ny)
    assert(z < nz)
    coordsToIdx(x+dx, y+dy, z+dz)
  }
  val xp = mkOffsetArray(+1, 0, 0)
  val xm = mkOffsetArray(-1, 0, 0)
  val yp = mkOffsetArray(0, +1, 0)
  val ym = mkOffsetArray(0, -1, 0)
  val zp = mkOffsetArray(0, 0, +1)
  val zm = mkOffsetArray(0, 0, -1)
  
  // central derivatives
  def d_x(f: Int => Double, i: Int) = (f(xp(i)) - f(xm(i))) / (2*dx)
  def d_y(f: Int => Double, i: Int) = (f(yp(i)) - f(ym(i))) / (2*dx)
  def d_z(f: Int => Double, i: Int) = (f(zp(i)) - f(zm(i))) / (2*dx)
  def laplacian(f: Int => Double, i: Int) =
    (  f(xp(i)) + f(xm(i))
     + f(yp(i)) + f(ym(i))
     + f(zp(i)) + f(zm(i))
     - 6 * f(i)) / (dx*dx)
  def laplacian2(f: Int => Double, i: Int) =
    (  f(xp(xp(i))) + f(xm(xm(i)))
     + f(yp(yp(i))) + f(ym(ym(i)))
     + f(zp(zp(i))) + f(zm(zm(i)))
     - 6 * f(i)) / (4*dx*dx)
  
  def hypot(x: Double, y: Double, z: Double) = math.sqrt(x*x + y*y + z*z)
  def sqr(x: Double) = x*x
  def cube(x: Double) = x*x*x
  
  // contract two symmetric tensors
  def dotdot(t1: symten, t2: symten) = {
    t1.xx*t2.xx + t1.yy*t2.yy + t1.zz*t2.zz + 2*(t1.xy*t2.xy + t1.yz*t2.yz + t1.zx*t2.zx) 
  }
  
  init()

  def init() {
    for (i <- 0 until N) {
      u(i).x = 0
      u(i).y = 0
      u(i).z = 0
      udot(i).x = 0
      udot(i).y = 0
      udot(i).z = 0
      
      // insert martensitic nucleus
      val (x, y, z) = idxToCoords(i)
      val r = dx * hypot(x-nx/2, y-ny/2, z-nz/2)
      for (p <- 0 until 3)
        eta(p)(i) = 0
      if (r < L_nuc) {
        eta(0)(i) = 1
      }
    
//      eta(0)(i) = 1
//      eta(1)(i) = 1
//      eta(2)(i) = 1
      
      // u(i).x = math.exp(-(x-nx/2)*(x-nx/2)/20.0)*dx
    }
    
    buildTransfStrain()
//    for (i <- 0 until N) {
//      eps_pl(i) = symten(-eps_tr(i).xx, -eps_tr(i).yy, -eps_tr(i).zz, -eps_tr(i).xy, -eps_tr(i).yz, -eps_tr(i).zx)
//    }
  }

  
  def buildTransfStrain() {
    for (i <- 0 until N) {
      eps_tr(i) = symten(0, 0, 0, 0, 0, 0)
      for (p <- 0 until 3) {
        eps_tr(i).xx += eta(p)(i) * bain(p).xx
        eps_tr(i).yy += eta(p)(i) * bain(p).yy
        eps_tr(i).zz += eta(p)(i) * bain(p).zz
        eps_tr(i).xy += eta(p)(i) * bain(p).xy
        eps_tr(i).yz += eta(p)(i) * bain(p).yz
        eps_tr(i).zx += eta(p)(i) * bain(p).zx
      }
    }
  }
  
  def step() {
    
    // build transformation strain from fields eta
    buildTransfStrain()
    
    // build total strain from displacement vector u using central differences
    for (i <- 0 until N) {
      eps_tot(i).xx = d_x(u(_).x, i)
      eps_tot(i).yy = d_y(u(_).y, i)
      eps_tot(i).zz = d_z(u(_).z, i)
      eps_tot(i).xy = 0.5 * (d_x(u(_).y, i) + d_y(u(_).x, i))
      eps_tot(i).yz = 0.5 * (d_y(u(_).z, i) + d_z(u(_).y, i))
      eps_tot(i).zx = 0.5 * (d_z(u(_).x, i) + d_x(u(_).z, i))
    }
    
    // build elastic strain from other strains
    for (i <- 0 until N) {
      eps_el(i).xx = eps_tot(i).xx - eps_tr(i).xx - eps_pl(i).xx
      eps_el(i).yy = eps_tot(i).yy - eps_tr(i).yy - eps_pl(i).yy
      eps_el(i).zz = eps_tot(i).zz - eps_tr(i).zz - eps_pl(i).zz
      eps_el(i).xy = eps_tot(i).xy - eps_tr(i).xy - eps_pl(i).xy
      eps_el(i).yz = eps_tot(i).yz - eps_tr(i).yz - eps_pl(i).yz
      eps_el(i).zx = eps_tot(i).zx - eps_tr(i).zx - eps_pl(i).zx
    }
    
    val c0 = E / (1+nu)
    val c1 = nu / (1 - 2*nu)
    
    // build sigma from elastic strain according to isotropic elasticity
    for (i <- 0 until N) {
      val eps_kk = eps_el(i).xx + eps_el(i).yy + eps_el(i).zz
      sigma(i).xx = c0 * (eps_el(i).xx + c1*eps_kk)
      sigma(i).yy = c0 * (eps_el(i).yy + c1*eps_kk)
      sigma(i).zz = c0 * (eps_el(i).zz + c1*eps_kk)
      sigma(i).xy = c0 * (eps_el(i).xy)
      sigma(i).yz = c0 * (eps_el(i).yz)
      sigma(i).zx = c0 * (eps_el(i).zx)
    }
    
    // calculate the divergence of the stress using central differences, but
    // with special treatment for second derivatives
    for (i <- 0 until N) {
      divsig(i).x = d_x(sigma(_).xx, i) + d_y(sigma(_).xy, i) + d_z(sigma(_).zx, i)
      divsig(i).y = d_x(sigma(_).xy, i) + d_y(sigma(_).yy, i) + d_z(sigma(_).yz, i)
      divsig(i).z = d_x(sigma(_).zx, i) + d_y(sigma(_).yz, i) + d_z(sigma(_).zz, i)
      
      // add and subtract different laplacian stencils for improved accuracy
      divsig(i).x += c0 * (1 + c1) * (laplacian(u(_).x, i) - laplacian2(u(_).x, i))
      divsig(i).y += c0 * (1 + c1) * (laplacian(u(_).y, i) - laplacian2(u(_).y, i))
      divsig(i).z += c0 * (1 + c1) * (laplacian(u(_).z, i) - laplacian2(u(_).z, i))
    }
    
    // integrate displacement vector one timestep using velocity verlet
    for (i <- 0 until N) {
      udot(i).x += dt * (divsig(i).x/rho - udot(i).x/gamma0) 
      udot(i).y += dt * (divsig(i).y/rho - udot(i).y/gamma0) 
      udot(i).z += dt * (divsig(i).z/rho - udot(i).z/gamma0) 
      
      u(i).x += dt * udot(i).x
      u(i).y += dt * udot(i).y
      u(i).z += dt * udot(i).z
    }
    
    // integrate relaxation dynamics of eta_p one timestep
    for (i <- 0 until N) {
      for (p <- 0 until 3) {
        laplace_eta(p)(i) = laplacian(eta(p)(_), i)
      }
    }
    for (i <- 0 until N) {
      val eta2 = sqr(eta(0)(i)) + sqr(eta(1)(i)) + sqr(eta(2)(i))
      for (p <- 0 until 3) {
        val etap = eta(p)(i)
        val dGgl_deta = -beta*laplace_eta(p)(i) + 2*A_hat*etap - 3*B_hat*sqr(etap) + 4*C_hat*eta2*etap
        val dGel_deta = (
           - sigma(i).xx * bain(p).xx
           - sigma(i).yy * bain(p).yy
           - sigma(i).zz * bain(p).zz
           - 2 * sigma(i).xy * bain(p).xy
           - 2 * sigma(i).yz * bain(p).yz
           - 2 * sigma(i).zx * bain(p).zx
        )
        eta(p)(i) -= dt * M * (dGgl_deta + dGel_deta)
      }
    }
    
    // integrate relaxation dynamics of plastic strain one timestep
//    for (i <- 0 until N) {
//      val e_kk = 0 * (eps_el(i).xx + eps_el(i).yy + eps_el(i).zz)
//      eps_pl(i).xx -= dt * k * (eps_el(i).xx - e_kk/3.0)  
//      eps_pl(i).yy -= dt * k * (eps_el(i).yy - e_kk/3.0)  
//      eps_pl(i).zz -= dt * k * (eps_el(i).zz - e_kk/3.0)  
//      eps_pl(i).xy -= dt * k * (eps_el(i).xy)  
//      eps_pl(i).yz -= dt * k * (eps_el(i).yz)  
//      eps_pl(i).zx -= dt * k * (eps_el(i).zx)
//    }
 
    // energy densities
    for (i <- 0 until N) {
      el_energy(i) = 0.5 * dotdot(sigma(i), eps_el(i))
      
      gl_energy(i) = 0.0
      val eta2 = sqr(eta(0)(i)) + sqr(eta(1)(i)) + sqr(eta(2)(i))
      for (p <- 0 until 3) {
        val etap = eta(p)(i)
        gl_energy(i) += -0.5*beta*etap*laplace_eta(p)(i) + A_hat*sqr(etap) - B_hat*cube(etap)
      }
      gl_energy(i) += C_hat*sqr(eta2)
    }
    timeSteps += 1 
  }
  
  def time = dt * timeSteps
  
  def energy() = gl_energy.sum + el_energy.sum
  
}

object YedduFccBctApp extends App {
  new Control(new YedduFccBctApp(), "Fcc to Bct (Kachaturyan)")  
}

class YedduFccBctApp extends Simulation {
  val gridEta0 = new Grid("eta0")
  val gridEta1 = new Grid("eta1")
  val gridEta2 = new Grid("eta2")
  val gridEl   = new Grid("elastic energy")
  val gridGl   = new Grid("ginzburg-landau energy")
  val plot    = new Plot("cross")
  
  var sim: FccBctSim = _
  
  def load(c: Control) {
    c.frameTogether("grids", gridEta0, gridEta1, gridEta2, gridEl, gridGl)
    
    params.add("L (um)", 1.0)
    params.add("L_nuc (um)", 0.1) 
    params.add("n0", 50)
    params.add("Dg (J/m^3)", -7e8)
    params.add("gamma (J/m^2)", 0.01)
    params.add("delta / dx", 2.0)
    params.add("E (GPa)", 200)
    params.add("nu", 0.3)
    params.add("b3", -0.1882)
    params.add("b1", 0.1363)
    params.add("rho (kg/m^3)", 7.8e3)
    params.add("dt * (v / dx)", 0.2)
    params.add("sigma_y (MPa)", 800.0)
    params.add("gamma0 / tau", 10.0)
    params.add("L delta Dg / v", 0.5)
    params.add("k / tau", 1.0)
    
    params.add("energy")
    params.add("time")
    
    params.addm("slice", new DoubleValue(0, 0, 1).withSlider)
  }
  
  
  def animate() {
    val sliceIdx = math.min((params.fget("slice") * sim.nz).toInt, sim.nz-1)
    def slice(a: Array[Double]) = {
      a.slice(sliceIdx*sim.nx*sim.ny, (sliceIdx+1)*sim.nx*sim.ny)
    } 
    
    val n0 = sim.n0
    gridEta0.registerData(n0, n0, slice(sim.eta(0)))
    gridEta1.registerData(n0, n0, slice(sim.eta(1)))
    gridEta2.registerData(n0, n0, slice(sim.eta(2)))
    gridEta0.setScale(0, 1)
    gridEta1.setScale(0, 1)
    gridEta2.setScale(0, 1)
    
    gridEl.registerData(n0, n0, slice(sim.el_energy))
    gridGl.registerData(n0, n0, slice(sim.gl_energy))

//    val ux = sim.u.map(_.x)
//    val ux_x = Array.tabulate(sim.N)(i => sim.d_x(ux(_), i))
//    val uy_y = Array.tabulate(sim.N)(i => sim.d_y(sim.u.map(_.y), i))
//    plot.registerLines("ux_x", new PointSet(0, sim.dx, ux_x.slice(0, sim.nx)), java.awt.Color.RED)
    
    params.set("energy", sim.energy)
    params.set("time", sim.time)
  }

  def clear() {
    gridEta0.clear()
    gridEta1.clear()
    gridEta2.clear()
    gridEl.clear()
    gridGl.clear()
    plot.clear()
  }

  def run() {
    sim = new FccBctSim(params)

    Job.animate();
    while (true) {
      for (i <- 0 until 1)
        sim.step();
      Job.animate();
      
//      if (sim.time > 2.5e-10)
//        Job.signalStop()
    }
  }
}
