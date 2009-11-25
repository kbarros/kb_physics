#!/usr/bin/python

import os
from itertools import *
from math import *


# -----------------------------------------------------------------
# Create a new directory, with a possible suffix for uniqueness.
# Return the name of the newly created directory.
#

def mknewdir(newdir):
    if os.path.exists(newdir):
        suffix = 1
        while os.path.exists(newdir+"-"+str(suffix)):
            suffix += 1
        newdir += "-"+str(suffix)
    os.mkdir(newdir)
    return newdir


# -----------------------------------------------------------------
# Simulation configuration
#

#
# NOTES:
# 1)  for Jusufi's parameters, energy of single chain with N=2 is
#   1.21375, where there is a single FENE bond, coloumbic attraction,
#   and double LJ repulsion, all at r=1
#   more details :
#    * total electrostatic energy is -1/\eps = -3.6  (per particle is -1.8)
#    * total double LJ repulsion energy is 2  (per particle is 1) --- LJ shift is +1
#
#   virial contribution to pressure is (1/3V) ( < r_ij F_ij > + < U_coul > )
#     where < r_ij F_ij > = LJ_force = 24   (note: globe forces don't contribute to virial!)
#           < U_coul > = -3.6
#     (note: sum runs over pairwise interactions with i<j)
#
# 2)  to test against pai-yi's data, the following changes are necessary:
#   ewald_cutoff=10, N=16, f=1, Ns=??, L=12.5992, salt_free=False, globe_epsilon=0
#

ewald_accuracy = 1e-4
ewald_cutoff = 20 # tunable parameter for performance
packing_distance = 1.0 # separation of ions when initially packed

equilibsteps = int(5e5)
simsteps     = int(1e7)
dumpevery = int(1e4)
num_dumps = simsteps / dumpevery

sim_hours = 5000

dt = 0.003 # time step // 0.002
N = 30 # chain length // 30
f = 40 # number of chains // 40
Ns = 400 # salt molecules
Z = 3 # salt valency // 1 through 4
R = 6 # globe radius // 6
L = 92 # linear system size // 160
T = 1.2 # temperature // 1.2
bjerrum = 3.0 # bjerrum length // 3.0

sigma_mm = 1.0 # monomer diameter for LJ repulsion // 1.0
sigma_ii = 1.0 # ion diameter for LJ repulsion // 1.0
sigma_mi = (sigma_mm + sigma_ii) / 2.0


salt_free = True # if True, only counterions (no salt coions)
if salt_free:
    n_monomers = f*N
    n_counterions = f*N - Z*Ns
    n_salt_counterions = Ns
    n_salt_coions = 0
    jobname = "L%d.Ns%d.sf.r9" % (round(L), Ns) # -round(log10(ewald_accuracy))
else:
    n_monomers = f*N
    n_counterions = f*N
    n_salt_counterions = Ns
    n_salt_coions = Z*Ns
    jobname = "L%d.Ns%04d.r8" % (round(L), Ns) # -round(log10(ewald_accuracy))

if n_counterions < 0:
    print "Too many multivalent ions"
    exit()
elif n_counterions == 0:
    print "Phi ratio: inf"
else:
    print "Phi ratio: %f" % (n_salt_counterions / float(n_counterions))

dirname = mknewdir("/home/kbarros/scratch/peb/"+jobname)
print "Created directory:\n " + dirname
lammps = "/home/kbarros/installs/Lammps/current/src/lmp_suse_linux"

confname = "in.dat"
dataname = "data.dat"
chaincfgname = "chaincfg.dat"
analyzername = "analyzer.scala"
cmdsname = "cmds.dat"
dumpname = "dump.dat"
pbsname = "job.pbs"



# ------------------------------------------------------------------
# This arrangement of points was acquired from [1]
# [1] R. H. Hardin, N. J. A. Sloane and W. D. Smith, Tables of putatively
# optimal packings on the sphere, published electronically at
# http://www.research.att.com/~njas/packings/
def load_sphere_packing():
    # combine list items into groups of size n
    def grouper(n, iterable, padvalue=None):
        "grouper(3, 'abcdefg', 'x') --> ('a','b','c'), ('d','e','f'), ('g','x','x')"
        return izip(*[chain(iterable, repeat(padvalue, n-1))]*n)
    
    with open("sphere_packings/pack.3.%d.txt"%f) as file:
        coords = map(float, file.readlines())
        return [tuple(e) for e in grouper(3, coords)]

sphere_packing = load_sphere_packing()


# -----------------------------------------------------------------
# Generate LAMMPS config file
#

def generate_conf_file():
    return (
"""# automatically generated polymer configuration file

units		lj
boundary	p p p
dimension	3

atom_style	full
neighbor	0.3 bin
neigh_modify	delay 5 page 500000 one 10000 # defaults: page 100000 one 2000

read_data	%(dataname)s		# load volume dimensions, masses, atoms, and bonds
# read_restart	restart.equil.dat
restart		%(restart_freq)d restart/restart.*.dat
thermo		%(dumpevery)d		# steps between output of thermodynamic quantities 
thermo_modify	flush yes		# flush output in case of Lammps crash

### BONDS
bond_style	fene
bond_coeff	* 7 2 0 0		# < bond_type K R0 fene_epsilon fene_sigma >
special_bonds	lj/coul 1.0 1.0 1.0  	# turn on regular interactions between all monomers (even bonded ones)
timestep	%(dt)f                  # FENE bonds are bottleneck for dt

### DEFINE GROUPS
group		graft type 1
group		nongraft subtract all graft
group		chain type 2
group		ion subtract all graft chain

### ADD GLOBE
fix		1 chain globe/lj126 %(R)s 1.0 %(sigma_mm)s %(cut_mm)s # radius epsilon sigma cutoff
fix_modify	1 energy yes		# include globe energy
fix		2 ion globe/lj126 %(R)s 1.0 %(sigma_mi)s %(cut_mi)s # radius epsilon sigma cutoff
fix_modify	2 energy yes		# include globe energy

### SET INTEGRATION METHOD
fix		3 nongraft nve		# apply newtons equations to ungrafted atoms
fix		4 graft nve/noforce	# zero velocities of grafted atoms while retaining force information

#### SOFT DYNAMICS TO SPREAD ATOMS
pair_style	soft 1.0			# < cutoff  >
pair_coeff	* * 0.0 60.0 			# < Astart Astop >
fix		5 nongraft langevin %(T)s %(T)s 10.0 699483	# < Tstart Tstop damp seed [keyword values ... ] >
run		%(soft_equilibsteps)d
unfix		5

### SWITCH TO COULOMB/LJ INTERACTIONS
pair_style	lj/cut/coul/long 0 %(ewald_cutoff)f	# < LJ_cutoff coulomb_cutoff > (LJ_cutoff is overridden below)
pair_coeff	1*2 1*2 1.0 %(sigma_mm)s %(cut_mm)s	# < type1 type2 epsilon sigma LJ_cutoff > (monomer-monomer)
pair_coeff	1*2 3*  1.0 %(sigma_mi)s %(cut_mi)s	# < type1 type2 epsilon sigma LJ_cutoff > (monomer-ion)
pair_coeff	3* 3*   1.0 %(sigma_ii)s %(cut_ii)s	# < type1 type2 epsilon sigma LJ_cutoff > (ion-ion)

pair_modify	shift yes		# LJ interactions shifted to zero
kspace_style	pppm %(ewald_accuracy)f	# desired accuracy
# LAMMPS coulombic energy is (q1 q2) / (\eps r)
# bjerrum length is (\lambda_B = %(bjerrum)f)
# temperature is (k T = %(T)f)
# dielectric is (1 / k T \lambda_B)
dielectric	%(dielectric)s

### EQUILIBRATE WITH TEMPERATURE RESCALING
fix		6 all temp/rescale 1 %(T)s %(T)s 0.05 1	# N Tstart Tstop window fraction
run		%(equilibsteps)d
unfix		6

### PRODUCTION RUN WITH LANGEVIN DYNAMICS
dump		1 all custom %(dumpevery)d %(dumpname)s id type x y z ix iy iz vx vy vz
fix		7 nongraft langevin %(T)s %(T)s 10.0 699483	# < Tstart Tstop damp seed [keyword values ... ] >
run		%(simsteps)d
unfix		7
""" % {'dataname':dataname, 'dumpname':dumpname, 'simsteps':simsteps,
       'soft_equilibsteps':equilibsteps, 'equilibsteps':equilibsteps, 'dumpevery':dumpevery, 'restart_freq':20*dumpevery,
       'ewald_accuracy':ewald_accuracy, 'ewald_cutoff':ewald_cutoff, 'T':T, 'R':R, 'dt':dt, 'bjerrum':bjerrum,
       'dielectric':1.0/(T * bjerrum),
       'sigma_mm':sigma_mm, 'sigma_mi':sigma_mi, 'sigma_ii':sigma_ii,
       'cut_mm':sigma_mm*(2**(1./6)), 'cut_mi':sigma_mi*(2**(1./6)), 'cut_ii':sigma_ii*(2**(1./6))})


# -----------------------------------------------------------------
# Generate input data file
#

def generate_data_file():
    natoms = n_monomers + n_counterions + n_salt_counterions + n_salt_coions
    
    header = (
"""LAMMPS FENE chain data file

          %(natoms)d  atoms
          %(nbonds)d  bonds
           0  angles
           0  dihedrals
           0  impropers

           5  atom types
           1  bond types
           0  angle types
           0  dihedral types
           0  improper types

 %(lmin)f  %(lmax)f xlo xhi
 %(lmin)f  %(lmax)f ylo yhi
 %(lmin)f  %(lmax)f zlo zhi

""" % {'natoms':natoms, 'nbonds':f*(N-1), 'lmin':0, 'lmax':L})
    
    masses = (
"""Masses
  #
  1  1.0 # grafted monomers
  2  1.0 # monomers
  3  1.0 # counterions
  4  1.0 # salt counterions
  5  1.0 # salt coions

""")
    
    def wrap_coordinate(x):
        i = 0
        while x >= L:
            x -= L
            i += 1
        while x < 0:
            x += L
            i -= 1
        return (x, i)
    
    def atom_position(i):
        d = packing_distance
        mx = d
        my = d*d / L
        mz = d*d*d / L**2
        x = i*mx
        y = i*my
        z = i*mz
        return (x%L,y%L,z%L)
    
    def polymer_position(chain, len):
        return [(len+R)*x+L/2. for x in sphere_packing[chain]]
    
    def atom_str((x,y,z), atom_id, molecule_id, atom_type, charge):
        (x,ix), (y,iy), (z,iz) = map(wrap_coordinate, (x,y,z))
        return "  %d  %d  %d  %f  %f  %f  %f  %d  %d  %d\n" % (atom_id, molecule_id, atom_type, charge, x, y, z, ix, iy, iz)
    
    atoms = "Atoms\n  # atom_id, molecule_id, atom_type, charge, x, y, z, wind_x, wind_y, wind_z \n"
    
    atom_id = 0
    
    # create polymer chain with f*N monomers
    for chain in range(f):
        for i in range(N):
            atom_id += 1
            molecule_id = chain+1
            atom_type = 1 if i == 0 else 2
            charge = -1
            atoms += atom_str(polymer_position(chain, i), atom_id, molecule_id, atom_type, charge)
    # create counter ions
    for i in range(n_counterions):
        atom_id += 1
        molecule_id = -1
        atom_type = 3
        charge = 1
        atoms += atom_str(atom_position(atom_id), atom_id, molecule_id, atom_type, charge)
    # create salt counterions
    for i in range(n_salt_counterions):
        atom_id += 1
        molecule_id = -1
        atom_type = 4
        charge = Z
        atoms += atom_str(atom_position(atom_id), atom_id, molecule_id, atom_type, charge)
    # create salt coions
    for i in range(n_salt_coions):
        atom_id += 1
        molecule_id = -1
        atom_type = 5
        charge = -1
        atoms += atom_str(atom_position(atom_id), atom_id, molecule_id, atom_type, charge)
    atoms += "\n"
    
    bonds = "Bonds\n  # bond_id, molecule_id, atom1_id, atom2_id\n"
    bond_id = 0
    for chain in range(f):
        for i in range(N-1):
            bond_id += 1
            bond_type = 1
            atom1 = N*chain + i+1
            atom2 = atom1 + 1
            bonds += "  %d  %d  %d  %d\n" % (bond_id, bond_type, atom1, atom2)
    bonds += "\n"
    
    return header + masses + atoms + bonds


# -----------------------------------------------------------------
# Generate configuration file for chain analysis
#

def generate_chain_config_file():
    return (
"""
dumpname       = %s
timesteps      = %d
maxdeltat      = 1
chainlength    = %d
no_chains      = %d
ad_condition_z = 1 
compute_msd    = 0 
outputfile     = analysis.dat
adsorb         = adsorb.dat 
msdout         = msd.dat
""" % (dumpname, num_dumps, N, f))


# -----------------------------------------------------------------
# Generate scala file for data analysis
#

def generate_scala_analyzer_file():
    return (
""" // run using command "sc analyzer.scala"

import kip._

val filename = "%s"
val chainLength = %d
val numChains = %d
val coreRadius = %f
PEB.go(chainLength, numChains, coreRadius)
""" % (dumpname, N, f, R))


# -----------------------------------------------------------------
# Generate commands file
#

def generate_cmds_file():
    return (
"""%(lammps)s < %(confname)s
sc analyzer.scala
generic_analyzer gyr.dat 3
""" % {'lammps':lammps, 'confname':confname})


# -----------------------------------------------------------------
# Generate PBS file
#

def generate_pbs_file():
    return (
"""
# ### AUTOMATICALLY GENERATED BATCH FILE

# ### name of job
#PBS -N %(jobname)s

# ### mail on end/abort
#PBS -m ea
#PBS -M kbarros@northwestern.edu

# ### maximum cpu time
#PBS -l cput=%(sim_hours)d:00:00

# ### number of nodes and processors per node
#PBS -l nodes=1:ppn=1

# ### indicates that job should not rerun if it fails
# #PBS -r n

# ### stdin and stderr merged as stderr
#PBS -j eo

# ### write stderr to file
#PBS -e %(dirname)s/log.err

# ### the shell that interprets the job script
#PBS -S /bin/bash

cd %(dirname)s 
%(lammps)s < %(confname)s > lammps.out
ssh minotaur scp -r %(dirname)s achilles.ms.northwestern.edu:%(dirname)s-mino
if [ $? -eq 0 ] ; then
touch COMPLETED
fi
""" % {'jobname':jobname, 'sim_hours':sim_hours, 'dirname':dirname, 'lammps':lammps, 'confname':confname})



# -----------------------------------------------------------------
# Build a new simulation directory
#

def build_dir():
    root = dirname + "/"
    
    os.mkdir(root+"restart")    # make directory for restarts
    
    with open(root + confname, 'w') as f:
        f.write(generate_conf_file())
        
    with open(root + dataname, 'w') as f:
        f.write(generate_data_file())
        
    with open(root + chaincfgname, 'w') as f:
        f.write(generate_chain_config_file())
        
    with open(root + analyzername, 'w') as f:
        f.write(generate_scala_analyzer_file())
        
    with open(root + cmdsname, 'w') as f:
        f.write(generate_cmds_file())
        
    with open(root + pbsname, 'w') as f:
        f.write(generate_pbs_file())
        
    os.system('cp ' + os.path.realpath(__file__) + ' ' + root + 'generation.py')
    os.system('ln -fsn %s link##' % root)

build_dir()
