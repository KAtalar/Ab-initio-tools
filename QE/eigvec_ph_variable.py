#!/bin/python

###########################################################
# Script to create the 'pw.x' input file corresponding
# to a structure distorted due to a soft phonon mode
# (A more general version)
###########################################################
# Written by Kemal Atalar (July 16, 2019)
###########################################################

from __future__ import division
import numpy as np
import os
import sys


####### INPUT PARAMETERS & DIRECTORIES ############
#element = 2 # 1 for Rb, 2 for K, 3 for Na
element = int(sys.argv[4])

#---- Cases for different materials ------
if element == 1:
  atom = 'Rb'; atom_lc = 'rb' # lowercase
elif element == 2:
  atom = 'K'; atom_lc = 'k'
elif element == 3:
  atom = 'Na'; atom_lc = 'na'
else:
  print 'ERROR: Use a valid element number input'
#-----------------------------------------

ph_dir = sys.argv[2] # Phonon directory, e.g. phonons_112
pref = sys.argv[1] # Prefix, e.g. 'rb20'
amp = float(sys.argv[3]) # Amplitude of distirtion in Angstroms

dyn_f = '../%s/%s.dyn1'%(ph_dir,pref)
eigv = 1 # eigenvector of the dynmat to use 

#atom_no = 20
atom_no = int(pref.split(atom_lc)[1])

debug=False
###################################################

# Read eigenvector from the dyn file
eigvec = [] # Soft mode eigvec of the dynamical matrix at G
count = 0
with open(dyn_f) as file:
  for line in file:
    if 'freq' in line.split():
      if count == 0: freq = float(line.split()[-2])
      count += 1
    elif count == eigv:
      splt = line.split()
      eigvec.append( (float(splt[1]), float(splt[3]), float(splt[5])) )

# Modify the SCF input according the eigvec distortion
scf_f = pref + '.scf.in'
scf_dist_f = pref + '.dist.scf.in' # SCF file to write distorted coord.
with open(scf_f,'r') as f:
  scf_lines = f.readlines()

count_p = 0
for ind,line in enumerate(scf_lines):
  splt = line.split()
  if (atom in splt) and (len(splt) == 4):
    pos = (float(splt[1]), float(splt[2]), float(splt[3]))
    new_pos = np.array(pos) + np.array(eigvec[count_p])*amp
    scf_lines[ind] = '%s       %.10f     %.10f     %.10f \n'%(atom,new_pos[0],new_pos[1],new_pos[2])
    count_p += 1

with open(scf_dist_f,'w+') as f:
  f.writelines(scf_lines)

# Checks
if debug:
  for line in scf_lines:
    if atom in line.split(): print line

if len(eigvec) != atom_no: print 'Error: Wrong eigvec length'
if freq > 0: print 'Warning: Not a soft mode'

