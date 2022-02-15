#!/bin/python

###########################################################
# Script to recalculate lambda_v values from the linewidth
# and frequency values (which is set to 0.00 by default
# for frequencies below certain value)
# This enables calculation of el-ph coupling for small
# frequency phonons (especially for phasons)
###########################################################
# Written by Kemal Atalar (Aug 16, 2019)
###########################################################

from __future__ import division
import numpy as np
import os
import glob
import fileinput

###########################################################
# Parameters
lmbd_mode = [1,2,3,4] # Modes, \nu to recalculate
		      # (especially for phason modes)
q_no = 1 # q-pt for the dynamical matrix file
print_b = False # Printing
###########################################################

# Read frequencies
freq = []
freq_GHz = []
if print_b: 
  print 'Dyn files'
  print glob.glob('*.dyn*')
for f in glob.glob('*.dyn%i'%q_no):
  i = int(f.split('dyn')[1])
  if print_b: print i
  pref = f.split('.')[0]
  scf_f = pref + '.scf.in'
  with open(f) as file:
    for line in file:
      if 'q' in line.split():
        qobj = line.split('(')[-1].split(')')[0]
        q_point = qobj.split()
      if 'freq' in line.split():
        w = float(line.split()[-2])
        freq.append(w) # in cm-1
        freq_GHz.append(w*29.9792458)
  if len(freq) != 0: 
    if print_b: print 'SYSTEM = %s'%pref
    if print_b: print '   '
  if freq:
    if print_b: print 'freq', freq
    command = 'freq%i = [q_point,freq_GHz]'%i
    exec(command)
    if print_b:
      prnt_command = 'print freq%i'%i
      exec(prnt_command)
  freq = []; freq_GHz = []

# Obtain details of system
#no_atoms = len(freq1[1])
no_qpoints = len(glob.glob('*.dyn1'))-1
if print_b: print no_atoms, no_qpoints

# Output format template
template = "     lambda({0:5d})= {1:7.4f}   gamma={2:8.2f} GHz"

# Read gamma and modify lambda
lambd = []
dos_iter = 0; lambd_iter = 0
for f in glob.glob('elph_dir/elph.inp_lambda.%i'%q_no):
  i = int(f.split('lambda.')[1])
  command_freq = 'freq_c = freq%i'%i
  exec(command_freq)
  for line in fileinput.input(f):
    if fileinput.filelineno() == 1: q_pt = line.split()[:3]
    if 'DOS' in line:
      if dos_iter > 0:
        command = 'lambd%i = [q_pt,lambd]'%i
        exec(command)
      lambd = []; lambd_iter = 0
      # Read the value of DOS at Ef
      nef = float(line.split('=')[1].split()[0])
      #print 'DOS= %f'%nef
    if 'lambda(' in line:
      gamma = float(line.split()[-2]) # linewidth in GHz
      ghz2ry = 3288463.807502993
      # Calculate the lambda value
      lambd_v = gamma*ghz2ry/(np.pi*nef*freq_c[1][lambd_iter]**2)
      if lambd_iter+1 in lmbd_mode:
	print template.format(lambd_iter+1,lambd_v,gamma)
      else:
	print line[:-1]
      lambd_iter += 1
    else:
      print line[:-1]
