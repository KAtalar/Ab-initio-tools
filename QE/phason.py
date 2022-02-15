#!/bin/python

###########################################################
# Script to determine phason modes from the 'ph.x' output
# (python phason.py)
###########################################################
# If standard deviation is zero along x&y but non-zero
# along z, it hints for a phason mode in the host-guest
# structure
###########################################################
# Written by Kemal Atalar (Aug 18, 2019)
###########################################################

from __future__ import division
import numpy as np
import os
import sys
import glob

###########################################################
####### INPUT PARAMETERS & DIRECTORIES ####################
###########################################################

element = 2 # 1 for Rb, 2 for K, 3 for Na

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

#pref = sys.argv[1] # Prefix, e.g. 'rb20'
pref = glob.glob('*.dyn1')[0].split('.')[0]
atom_no = int(pref.split(atom_lc)[1])
dyn_f = '%s.dyn1'%(pref)
dyn_asr_f = 'dynmat.out'
eigv = 6 # eigenvector of the dynmat to use 

debug=False
##########################################################

# Read eigenvector from the dyn file
eigvec = []
eig_asr = []

# Read ASR applied eigenvectors
count = 0
with open(dyn_f) as f:
  eig_sub = []
  for line in f:
    if 'freq' in line.split():
      if count != 0 and count <= eigv: eigvec.append(np.array(eig_sub))
      eig_sub = []
      count += 1
    elif count > 0 and count <= atom_no:
      splt = line.split()
      eig_sub.append( [float(splt[1]), float(splt[3]), float(splt[5])] )

if debug: print ' ----- First eigvec: \n', eigvec[0]
if debug: print ' ----- First eigvec (ASR): \n', eig_asr[0]

# Look at the eigvectors and try to identify the phason
for ind, vec in enumerate(eigvec):

  if debug:
    if ind==0: print vec

  print ' ---- Eigvec(%i)'%(ind+1)

  # Method 1: Standard deviation and mean of eigvec
  av_xyz = np.mean(vec,axis=0) # averages over columns
  std_xyz = np.std(vec,axis=0)
  print 'Average std along axes:  %.5f  %.5f  %.5f '%(std_xyz[0],std_xyz[1],std_xyz[2])
'''
  # Method 2: Comparison to ASR aplied eigvec
  for ind_a, vec_asr in enumerate(eig_asr):
    diff = np.sum((vec-vec_asr)**2.0)/len(vec)
    print ' * Diff to eig_asr(%i) is   %.5f'%(ind_a+1, diff)
'''
