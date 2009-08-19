#!/usr/bin/python

N = 8
L = 12.5992104989487 # 2000^(1/3)

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


atoms = "Atoms\n  # atom_id, molecule_id, atom_type, charge, x, y, z, wind_x, wind_y, wind_z \n"
# create polymer chain
for i in range(N):
    atom_id = i+1
    molecule_id = 1
    atom_type = 1
    charge = -1
    x = atom_id * 1.0
    y = L/2.0
    z = L/2.0
    # valid for 'full' type atoms
    atoms += ("  %d  %d  %d  %f  %f  %f  %f  0  0  0\n" % (atom_id, molecule_id, atom_type, charge, x, y, z))
# create ions
for i in range(N, 2*N):
    atom_id = i+1
    molecule_id = 2
    atom_type = 2
    charge = 1
    x = (atom_id-N) * 1.0
    y = L/2.0 + 1.0
    z = L/2.0
    # valid for 'full' type atoms
    atoms += ("  %d  %d  %d  %f  %f  %f  %f  0  0  0\n" % (atom_id, molecule_id, atom_type, charge, x, y, z))
atoms += "\n"


bonds = "Bonds\n  # bond_id, molecule_id, atom1_id, atom2_id\n"
for i in range(N-1):
    bond_id = i+1
    molecule_id = 1
    bonds += "  %d  %d  %d  %d\n" % (bond_id, molecule_id, bond_id, bond_id+1)
bonds += "\n"


print (header + masses + atoms + bonds)

