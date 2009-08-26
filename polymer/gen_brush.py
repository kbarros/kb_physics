#!/usr/bin/python

import os
from itertools import *


N = 4 # chain length // 30
f = 4 # number of chains // 40
Ns = 1 # salt molecules // ??
Z = 3 # salt valency // 1 through 4
R = 6 # globe radius // 6
L = 2.5*(R+N) # linear system size // 160

equilibsteps = 0 # 100000
simsteps     = 10000 # 10000000
dumpevery = 50
num_dumps = simsteps / dumpevery

dirname = "/home/kbarros/scratch/peb/N%d.Ns%d.Z%d.L%d" % (N, Ns, Z, L)
confname = "in.dat"
dataname = "data.dat"
chaincfgname = "chaincfg.dat"
cmdsname = "cmds.dat"
dumpname = "dump.dat"
analysisname = "analysis.dat"


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
neigh_modify	delay 5

read_data	%(dataname)s		# load volume dimensions, masses, atoms, and bonds

pair_style	lj/cut/coul/long 1.12246204830937 12 # LJ_cutoff=2^{1/6}  [ coulomb_cutoff ]
pair_coeff	* * 1.0 1.0		# < atom_type1 atom_type2 epsilon sigma [ LJ_cutoff coulomb_cutoff ] >
pair_modify	shift yes		# LJ interactions shifted to zero
kspace_style	ewald 1.0e-4		# desired accuracy
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

# apply newtons equations
fix		1 all nve
fix		2 all langevin 1.2 1.2 10.0 699483 # < ID group-ID langevin Tstart Tstop damp seed [keyword values ... ] >
fix		3 graft setforce 0.0 0.0 0.0	   # graft atoms don't move
fix		4 nongraft globe/lj126 %(R)f 1.0 1.0 1.12246204830937 # radius epsilon sigma cutoff

# run
timestep	0.003
thermo		%(dumpevery)d				# output thermodynamic quantities every 1000 steps

run		%(equilibsteps)d

dump		1 all atom %(dumpevery)d %(dumpname)s	# < ID group-ID style every_N_timesteps file args >
dump_modify	1 image yes
#dump_modify	1 image yes scale no

dump		2 all custom %(dumpevery)d veldump.dat vx vy vz

run		%(simsteps)d
""" % {'dataname':dataname, 'dumpname':dumpname, 'simsteps':simsteps,
       'equilibsteps':equilibsteps, 'dumpevery':dumpevery, 'R':R})



# -----------------------------------------------------------------
# Generate input data file
#

def generate_data_file():
    natoms = f*N + f*N + Ns + Z*Ns # monomers, counterions, salt counterions, salt coions
     
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
  3  1.0 # ions
  4  1.0 # salt counterions
  5  1.0 # salt coions

""")
    
    def atom_position(i):
        mx = 1
        my = (1 / L)
        mz = (1 / L**2)
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
    
    # create polymer chain
    for chain in range(f):
        for i in range(N):
            atom_id += 1
            molecule_id = chain+1
            atom_type = 1 if i == 0 else 2
            charge = -1
            atoms += atom_str(polymer_position(chain, i), atom_id, molecule_id, atom_type, charge)
    # create ions
    for i in range(f*N):
        atom_id += 1
        molecule_id = -1
        atom_type = 2
        charge = 1
        atoms += atom_str(atom_position(atom_id), atom_id, molecule_id, atom_type, charge)
    # create salt counterions
    for i in range(Ns):
        atom_id += 1
        molecule_id = -1
        atom_type = 3
        charge = Z
        atoms += atom_str(atom_position(atom_id), atom_id, molecule_id, atom_type, charge)
    # create salt coions
    for i in range(Z*Ns):
        atom_id += 1
        molecule_id = -1
        atom_type = 4
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
# Build a new simulation directory
#

def build_sim_dir():
    root = mknewdir(dirname) + "/"
    print "creating directory:\npushd " + root
    
    with open(root + confname, 'w') as f:
        f.write(generate_conf_file())
    
    with open(root + dataname, 'w') as f:
        f.write(generate_data_file())
    
    with open(root + chaincfgname, 'w') as f:
        f.write(generate_chain_config_file())

    with open(root + cmdsname, 'w') as f:
        f.write(generate_cmds_file())

build_sim_dir()
