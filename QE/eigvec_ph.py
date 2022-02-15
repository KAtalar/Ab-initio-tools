#!/bin/python

###########################################################
# Script to create the 'pw.x' input file corresponding
# to a structure distorted due to a soft phonon mode
###########################################################
# Written by Kemal Atalar (May 30, 2019)
###########################################################

from __future__ import division
import numpy as np
import os
import sys


####### INPUT PARAMETERS & DIRECTORIES ############

#ph_dir = 'phonons_112' # Phonon directory
ph_dir = sys.argv[2]
#pref = 'rb20' # Prefix
pref = sys.argv[1]
atom = 'Rb'
dyn_f = '../%s/%s.dyn1'%(ph_dir,pref)
#amp = 0.01 # Amplitude of distirtion in Angstroms
amp = float(sys.argv[3])

###########################################################

# Read eigenvector from the dyn file
eigvec = [] # Soft mode eigvec of the dynamical matrix at G
count = 0
with open(dyn_f) as file:
  for line in file:
    if 'freq' in line.split():
      if count == 0: freq = float(line.split()[-2])
      count += 1
    elif count == 1:
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
  if atom in line[0:2]:
    pos = (float(splt[1]), float(splt[2]), float(splt[3]))
    new_pos = np.array(pos) + np.array(eigvec[count_p])*amp
    scf_lines[ind] = '%s       %.10f     %.10f     %.10f \n'%(atom,new_pos[0],new_pos[1],new_pos[2])
    count_p += 1

#for line in scf_lines:
#  if atom in line.split(): print line

with open(scf_dist_f,'w+') as f:
  f.writelines(scf_lines)

# Checks
atom_no = float(pref[2:])
print atom_no
if len(eigvec) != atom_no: print 'Error: Wrong eigvec length'
if freq > 0: print 'Warning: Not a soft mode'
