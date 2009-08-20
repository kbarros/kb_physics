#!/usr/bin/python

N = 16
rho = 0.008
L = (N / rho)**(1./3.) # approx 12.6


header = ("""LAMMPS FENE chain data file
 
          %(natoms)d  atoms
          %(nbonds)d  bonds
           0  angles
           0  dihedrals
           0  impropers
 
           2  atom types
           1  bond types
           0  angle types
           0  dihedral types
           0  improper types

 %(lmin)f  %(lmax)f xlo xhi
 %(lmin)f  %(lmax)f ylo yhi
 %(lmin)f  %(lmax)f zlo zhi

""" % {'natoms':2*N, 'nbonds':(N-1), 'lmin':0, 'lmax':L})


masses = (
"""Masses
  #
  1  1.0
  2  1.0

""")



def wrap_coordinate(x):
    return (x % L, x / L) # todo: handle negative winding numbers

def atom_position(i):
    mx = 1
    my = (1 / L)
    mz = (1 / L**2)
    x = L/2. + i*mx
    y = L/2. + i*my
    z = L/2. + i*mz
    return (x,y,z)


def atom_str(atom_id, molecule_id, atom_type, charge):
    (x,ix), (y,iy), (z,iz) = map(wrap_coordinate, atom_position(atom_id))
    return "  %d  %d  %d  %f  %f  %f  %f  %d  %d  %d\n" % (atom_id, molecule_id, atom_type, charge, x, y, z, ix, iy, iz)


atoms = "Atoms\n  # atom_id, molecule_id, atom_type, charge, x, y, z, wind_x, wind_y, wind_z \n"
# create polymer chain
for i in range(N):
    atom_id = i+1
    molecule_id = 1
    atom_type = 1
    charge = -1
    atoms += atom_str(atom_id, molecule_id, atom_type, charge)
# create ions
for i in range(N):
    atom_id = N+(i+1)
    molecule_id = 2
    atom_type = 2
    charge = 1
    atoms += atom_str(atom_id, molecule_id, atom_type, charge)
atoms += "\n"


bonds = "Bonds\n  # bond_id, molecule_id, atom1_id, atom2_id\n"
for i in range(N-1):
    bond_id = i+1
    molecule_id = 1
    bonds += "  %d  %d  %d  %d\n" % (bond_id, molecule_id, bond_id, bond_id+1)
bonds += "\n"


print (header + masses + atoms + bonds)

