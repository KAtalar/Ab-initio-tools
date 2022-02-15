#!/bin/python

###########################################################
# Script to extract frequencies from 'ph.x' dynamical 
# matrix output, and print to standard output
###########################################################
# Written by Kemal Atalar (Aug 18, 2019)
###########################################################

from __future__ import division
import numpy as np
import os
import glob

# Choose files including frequency info
asr = int(raw_input('ASR corrected/not for dynamical matrix? (1/0)'))

if asr == 1:
  dyn_file = 'dynmat.out'
elif asr == 0:
  dyn_file = '*.dyn*'
else:
  print 'ERROR: Wrong input, put 1 or 0'

# Get frequencies
freq = []
for f in glob.glob(dyn_file):
  pref = f.split('.')[0]
  scf_f = pref + '.scf.in'
  with open(f) as file:
    for line in file:
      if 'q' in line.split():
        q_point = line
      if 'freq' in line.split():
        freq.append(float(line.split()[-2]))
  if len(freq) != 0: 
    print 'SYSTEM = %s'%pref
    print q_point
    for ind,w in enumerate(freq):
      print 'freq (%i) = %f [cm-1], %e [Hartree]'%(ind+1,w,w*4.5563e-6)
    print '   '
  freq = []
