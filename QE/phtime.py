#!/bin/python

###########################################################
# Script to read the details of time taken for a 'ph.x' 
# calculation, separated for different q-points
# (python phtime.py)
###########################################################
# Written by Kemal Atalar (Aug 6, 2019)
###########################################################

from __future__ import division
import numpy as np
import os
import glob

###########################################################
# Find the phonon or electron-phonon output file in the current directory
if not glob.glob('*.ph*out*'):
  file_n =  glob.glob('*.elph*out*')[-1]
else:
  file_n =  glob.glob('*.ph*out*')[0]
print file_n

scf_f = glob.glob('*scf.in')[0]
pref = scf_f.split('.')[0]

rep_no = -1; core_no = -1,
nspace = -1; npool = -1
switch = 1
###########################################################

# Read k-grid points from scf input
with open(scf_f) as ff:
  k_ind = -1
  for ind,line in enumerate(ff):
    if 'K_POINTS' in line.split():
      k_ind = ind+1
    if ind == k_ind:
      kgrid = line.split()[:3]

# Read info and times from ph output
with open(file_n) as f:
  count_sec = 0; count_iter = 0
  total_sec = 0; total_iter = 0
  count_q = -1
  iter_bef = 0

  for line in f:
    # Read data
    if 'irreducible' in line.split():
      rep_no = int(line.split('irreducible')[0].split()[-1])
    if 'processors' in line.split():
      core_no = int(line.split('processors')[0].split()[-1])
    if 'space' in line.split():
      nspace = int(line.split()[-1])
    if 'K-points' in line.split():
      npool = int(line.split()[-1])

    # Initial representations
    if 'PHONON' == line[5:11] and switch == 1:
      switch = 0
      # CPU time
      rep_cpu = line.split(':')[1].split('CPU')[0]
      # Wall time
      rep_wall = line.split('CPU')[1].split('WALL')[0]

      print '    '
      print 'INFO ---- System: %s'%(pref)
      print '     -- scf k-grid = %sx%sx%s'%(kgrid[0],kgrid[1],kgrid[2])
      print '     -- no. of irreps = %i'%(rep_no)
      print '     -- no. of cores = %i'%(core_no)
      print '     -- no. of pools (k-point paral.) = %i'%(npool)
      print '     -- no. of cores for space paral. = %i'%(nspace)

      print 'TIME'
      print ' -- Time taken for deciding irreps ='
      print '%s CPU, %s WALL'%(rep_cpu, rep_wall)

    # Iterations
    if ('secs' in line.split()) and ('iter' in line.split()):
      count_sec += 1

      # Time taken per iteration
      secs = float(line.split(':')[1].split()[0])
      if count_sec != 1:
        sec_diff = secs - sec_bef
        total_sec += sec_diff
      elif count_q > 0:
	sec_diff = secs - sec_bef
        print ' -- Time for precalculations (q%i) = %f secs'%(count_q, sec_diff)
      else: print ' -- Time for precalculations (q%i) = %f secs'%(count_q, secs)
      sec_bef = secs

      # Iterations
      iter = int(line.split()[2])
      if iter < iter_bef:
        count_iter += 1
        total_iter += iter_bef
	'''
        if count_iter == rep_no:
          it_per_rep = total_iter/count_iter
          print ' -- Iterations per rep (q%i) = %f'%(count_q,it_per_rep)
          count_iter = 0; total_iter = 0
	'''
      iter_bef = iter

    # Calculate and reset when a q point is finished
    if 'Computing' in line.split():
      if count_q > -1: # After the first q calculation ends
	# Check if correct no of irreps are calculated
        if count_iter != rep_no-1:
	  print 'Error: no match for rep_no (%i != %i)'%(count_iter,rep_no-1)

        # Iterations per irrep
        it_per_rep = (total_iter+iter_bef)/(count_iter+1)
        print ' -- Iterations per rep (q%i) = %f'%(count_q,it_per_rep)

	# Time taken per iteration
        sec_per_iter = total_sec/count_sec
	print ' -- Seconds per iteration (q%i) = %f secs'%(count_q,sec_per_iter)
        tot_time = sec_per_iter*it_per_rep*rep_no
	print ' -- Estimated total time for q%i = %f secs'%(count_q,tot_time)
	print '  (without initial irrep time)   = %f min'%(tot_time/60)
	print '                                 = %f hr'%(tot_time/3600)
        print '    '

      # Reset
      iter_bef = 0 # To prevent counting of this iter into next q point
      count_iter = 0; total_iter = 0
      count_sec = 0; total_sec = 0
      count_q += 1

  # Calculate the average for next q point if the calculation is not finished yet
  if count_iter < rep_no and count_iter != 0:
    it_per_rep = total_iter/count_iter
    print ' -- Iterations per rep (q%i) upto %i reps = %f'%(count_q,count_iter,it_per_rep)
    sec_per_iter = total_sec/count_sec
    print ' -- Seconds per iteration (q%i) = %f secs'%(count_q,sec_per_iter)
    tot_time = sec_per_iter*it_per_rep*rep_no
    print ' -- Estimated total time for q%i = %f secs'%(count_q,tot_time)
    print '  (without initial irrep time)   = %f min'%(tot_time/60)
    print '                                 = %f hr'%(tot_time/3600)
    print '    '
