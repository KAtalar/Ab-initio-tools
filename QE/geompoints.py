#!/bin/python

###########################################################
# Script to extract positions after a geometry relaxation
# calculation by "pw.x" from the output
# (python geompoints.py >> SCFFILE)
###########################################################
# Written by Kemal Atalar (May 30, 2019)
###########################################################

import numpy as np
import os

cwd = os.getcwd()
#appr = 78
appr = int(cwd.split('_atom')[0].split('/')[-1])

#file_n = '%i_atom/geom_ideal_stricter/rb%i.geom.out'%(appr, appr)
file_n = 'rb%i.geom.out'%(appr)

pos=[]; latt=[]
with open(file_n) as f:
    count=0
    # Get atom positions
    for ind,line in enumerate(f):
        if 'Rb' in line[0:2]:
            if 'Rb' == line.split()[0]:
                pos.append(line.split()[1:])
        # Get latt parameters
        if 'CELL_PARAMETERS' in line.split():
            alat = line.split()[-1].split(')')[0]
            count+=1
        else: 
            if count > 0:
                latt.append(line.split())
                count += 1
        if count == 4:
            count = 0

alat_bohr = float(alat)
alat = float(alat)*0.529
latt_par = np.array(latt[-3:]).astype(float)
pos_A = np.array(pos[-appr:]).astype(float)

pos_scaled = np.zeros(np.shape(pos_A))
for ind_r,row in enumerate(pos_A):
    for ind_c,col in enumerate(row):
        pos_scaled[ind_r,ind_c] = col/(alat*latt_par[ind_c,ind_c])

print 'alat=%f'%alat
print 'Cell_parameters (Angstrom)'
for i in latt_par:
    print i[0]*alat,i[1]*alat,i[2]*alat
print 'Cell_parameters (a.u.)'
for i in latt_par:
    print i[0]*alat_bohr,i[1]*alat_bohr,i[2]*alat_bohr
print 'Atom positions (scaled)'
for i in pos_scaled:
    print i[0],i[1],i[2]
print 'Atom positions (Angstrom)'
for i in pos_A:
    print 'Rb', i[0],i[1],i[2]
print 'ibrav= 6, celldm(1) =%f, celldm(3) = %f'%(alat_bohr*latt_par[0][0],latt_par[2][2]/latt_par[1][1])
