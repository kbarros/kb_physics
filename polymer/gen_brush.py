#!/usr/bin/python


import os
from itertools import *


# -----------------------------------------------------------------
# Simulation configuration
#

# TODO: don't hardcode Bjerrum length

dt = 0.003
ewald_accuracy = 1e-5

N = 30 # chain length // 30
f = 40 # number of chains // 40
Ns = 240 # salt molecules
Z = 3 # salt valency // 1 through 4
R = 6 # globe radius // 6
L = 160 # linear system size // 160
T = 1.2 # temperature // 1.2

equilibsteps = 100000
simsteps     = 1000000
dumpevery = 5000
num_dumps = simsteps / dumpevery

salt_free = True # if True, only counterions (no salt coions)
if not salt_free:
    n_monomers = f*N
    n_counterions = f*N
    n_salt_counterions = Ns
    n_salt_coions = Z*Ns
else:
    n_monomers = f*N
    n_counterions = f*N - Z*Ns
    n_salt_counterions = Ns
    n_salt_coions = 0

dirname = "/home/kbarros/scratch/peb/N%d.f%d.Ns%d.Z%d.L%d" % (N, f, Ns, Z, L)
if salt_free:
    print "Phi ratio: %f" % (n_salt_counterions / float(n_counterions))
    dirname += ".sf"
lammps = "/home/kbarros/Lammps/current/src/lmp_suse_linux"

confname = "in.dat"
dataname = "data.dat"
chaincfgname = "chaincfg.dat"
cmdsname = "cmds.dat"
dumpname = "dump.dat"
analysisname = "analysis.dat"
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

pair_style	lj/cut/coul/long 1.12246204830937 25 # LJ_cutoff=2^{1/6}  [ coulomb_cutoff ]
pair_coeff	* * 1.0 1.0		# < atom_type1 atom_type2 epsilon sigma [ LJ_cutoff coulomb_cutoff ] >
pair_modify	shift yes		# LJ interactions shifted to zero
kspace_style	pppm %(ewald_accuracy)f	# desired accuracy
# LAMMPS coulombic energy is (q1 q2) / (\eps r)
# bjerrum length is (\lambda_B = 3)
# temperature is (k T = 1.2)
# so dielectric is (1 / (1.2 * 3))
dielectric	0.2777777777777

bond_style	fene
bond_coeff	* 7 2 0 0		# < bond_type K R0 fene_epsilon fene_sigma >
special_bonds	lj/coul 1.0 1.0 1.0  	# turn on regular interactions between all monomers (even bonded ones)

group		graft type 1
group		nongraft subtract all graft

fix		1 nongraft nve				# apply newtons equations
fix		2 nongraft globe/lj126 %(R)f 1.0 1.0 1.12246204830937 # radius epsilon sigma cutoff
fix_modify	2 energy yes				# include globe energy

timestep	%(dt)f
thermo		%(dumpevery)d				# steps between output of thermodynamic quantities 
restart		%(restart_freq)d restart.*.dat

# dump data in format for visualization by xmovie
# dump		1 all atom 200 movie.dat		# < ID group-ID style every_N_timesteps file args >

# equilibrate the simulation with a temp/rescale for stability
fix		3 all temp/rescale 1 %(T)f %(T)f 0.05 1	# N Tstart Tstop window fraction
run		%(equilibsteps)d
unfix		3

dump		2 all atom %(dumpevery)d %(dumpname)s	# < ID group-ID style every_N_timesteps file args >
dump_modify	2 image yes scale no

# run the simulation with langevin dynamics for performance & robustness
fix		4 nongraft langevin %(T)f %(T)f 10.0 699483	# < ID group-ID langevin Tstart Tstop damp seed [keyword values ... ] >
run		%(simsteps)d
""" % {'dataname':dataname, 'dumpname':dumpname, 'simsteps':simsteps,
       'equilibsteps':equilibsteps, 'dumpevery':dumpevery, 'restart_freq':20*dumpevery,
       'ewald_accuracy':ewald_accuracy, 'T':T, 'R':R, 'dt':dt})



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
    
    def atom_position(i):
        d = 8.0
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
        return "  %d  %d  %d  %f  %f  %f  %f  0  0  0\n" % (atom_id, molecule_id, atom_type, charge, x, y, z)
    
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
"""dumpname       = %s
timesteps      = %d
maxdeltat      = 1
chainlength    = %d
no_chains      = %d
outputfile     = %s
""" % (dumpname, num_dumps, N, f, analysisname))


# -----------------------------------------------------------------
# Generate commands file
#

def generate_cmds_file():
    return (
"""
lamps < %s
chain_analysis chaincfg.dat
generic_analyzer analysis.dat 8
""" % confname)


# -----------------------------------------------------------------
# Generate PBS file
#

def generate_pbs_file():
    return (
"""
# ### AUTOMATICALLY GENERATED BATCH FILE

# ### name of job
#PBS -N test240

# ### mail for begin/end/abort
#PBS -m bea
#PBS -M kbarros@northwestern.edu

# ### direct stderr and stdout to files
#PBS -e jobtest.err
#PBS -o jobtest.log

# ### maximum cpu time
#PBS -l cput=100:00:00

# ### number of nodes and processors per node
#PBS -l nodes=1:ppn=1

# ### indicates that job should not rerun if it fails
# #PBS -r n

# ### the shell that interprets the job script
#PBS -S /bin/bash

hostname

echo 'hello world'

cd %(dirname)s 
%(lammps)s < %(confname)s > job.log 2> job.err
""" % {'dirname':dirname, 'lammps':lammps, 'confname':confname})



# -----------------------------------------------------------------
# Build a new simulation directory
#

def build_sim_dir():
    root = mknewdir(dirname) + "/"
    print "creating directory:\n " + root
    
    with open(root + confname, 'w') as f:
        f.write(generate_conf_file())
    
    with open(root + dataname, 'w') as f:
        f.write(generate_data_file())
    
    with open(root + chaincfgname, 'w') as f:
        f.write(generate_chain_config_file())

    with open(root + cmdsname, 'w') as f:
        f.write(generate_cmds_file())
    
    with open(root + pbsname, 'w') as f:
        f.write(generate_pbs_file())
    
    os.system('ln -fsn %s link##' % root)

build_sim_dir()
