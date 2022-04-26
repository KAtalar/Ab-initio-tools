#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 15:40:41 2022

Script to convert Stillinger-Weber potential parameters defined for GULP into
LAMMPS format; and create the potential file for LAMMPS calculation

The equations describing 2-body and 3-body interactions of the 
Stillinger-Weber potential are as followed:
    
    ############## GULP ###############
    Parameters: A, rho, B, rmax & K, theta0, rho12, rho13, rmax12, rmax13
    
    V2(rij) = A*( B/(rij**4.) - 1 ) * exp(rho / (rij-rmax))
    
    V3(theta_ijk) = K*exp( rho12 / (rij-rmax12) + rho13 / (rik-rmax13) ) 
                    * ( cos(theta_ijk) - costheta0 )**2.
    
    ############## LAMMPS ###############
    Parameters: epsilon, A_L, B_L, sigma, p, q, a & lambda, gamma, costheta0
    
    V2(rij) = epsilon*A_L*( B_L*sigma**p/(rij**p) - sigma**q/(rij**p) )
              * exp(sigma / (rij - a*sigma))
    
    V3(theta_ijk) = epsilon * lambda 
                    * exp( gamma*sigma / (rij - a*sigma) 
                           + gamma*sigma / (rik - a*sigma) ) 
                    * ( cos(theta_ijk) - costheta0 )**2.

- Assumes a theta cutoff for the three-body interaction in LAMMPS as 
described in the appendix of Jiang et al.[1]

[1] Y. Zhou and J. Jiang, Scientific Reports 7, (2017).

@author: K. Atalar
"""

import numpy as np
from datetime import date
import itertools
import sys

###############################################################################
# Global variables
PFORMAT = '{:<3} {:<3} {:<3} {:>8.3f} {:>7.3f} {:>7.3f} {:>7.3f} {:>7.3f} {:>7.3f} {:>7.3f} {:>7.3f} {:>2d} {:>2d} {:>3.1f}'
METAL= ['Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Zr','Nb','Mo','Tc','Rh', \
         'Pd','Sn','Hf','Ta','W','Re','Ir','Pt']
CHALCOGEN = ['O','S','Se','Te']
###############################################################################

def readGulpInput(fname):
    '''
    Function to read GULP parameters from the input file
    '''
    # Initialize the lists
    atomName = []
    gulp2 = []
    gulp3_mxx = []; gulp3_xmm = []; gulp3_mx1x2 = []
    twobody, threebody = False, False
    
    with open(fname) as f:
        for line in f:
            splt = line.split()
            
            if twobody and line != '\n':
                # Read atom types
                at1, at2 = splt[:2]
                if at1 not in atomName:
                    atomName.append(at1)
                if at2 not in atomName:
                    atomName.append(at2)
                    
                # Read 2body parameters
                if len(gulp2) == 0:
                    gulp2 = splt[2:5]+[splt[6]]
                    gulp2 = [float(i) for i in gulp2]

                    
            elif threebody and line != '\n':
                gulp3 = splt[3:7]+[splt[8],splt[10]]
                gulp3 = [float(i) for i in gulp3]
                at3 = splt[:3]
                if at3[0] in CHALCOGEN or at3[0][:-1] in CHALCOGEN:
                    if len(gulp3_xmm) == 0:
                        gulp3_xmm = gulp3[:]
                elif at3[0] in METAL:
                    if at3[1] == at3[2]:
                        if len(gulp3_mxx) == 0:
                            gulp3_mxx = gulp3[:]
                    else:
                        if len(gulp3_mx1x2) == 0:
                            gulp3_mx1x2 = gulp3[:]
                            
                else:
                    print('Error in readGulpInput: Cannot identify whether atom is metal or chalcogen')
                
            
            # Control the switches for reading the 2body and 3body parameters
            if line != '\n' and splt[0] == 'sw2': 
                twobody=True
            elif line != '\n' and splt[0] == 'sw3':
                threebody=True
            elif line == '\n' and (threebody or twobody):
                twobody, threebody = False, False

    # Return the parameters
    if len(atomName) == 3:
        return atomName, [gulp2, gulp3_mxx, gulp3_xmm, gulp3_mx1x2]
    
    elif len(atomName) == 2:
        return atomName, [gulp2, gulp3_mxx, gulp3_xmm]
    
    else:
        print("Error in readGulpInput: Unexpected number of atoms - check GULP input file")
        sys.exit()
        
       
def gulpDic(atnam, gulp2, gulp3mxx, gulp3xmm, gulp3mx1x2=[], gulp2mm=[], gulp2xx=[], partype = 'Jiang'):
    '''
    Function to calculate the possible variations of interactions
    for a single structure and saving it to a dictionary
    '''
    # Define a zero gulp3 for Kandemir 2-body parameters
    gulp3zero = gulp3mxx.copy()
    gulp3zero[0] = 0.0 # Set K to zero for the 3-body interaction
    
    outdic = {}
    if len(atnam) == 2:
        mm, xx = atnam
        outdic[''.join([mm,xx,xx])] = gulp2 + gulp3mxx
        outdic[''.join([xx,mm,mm])] = gulp2 + gulp3xmm
        # Add other 2body interactions if Kandemir
        if partype[0] == 'K' or partype[0] == 'k':
            outdic[''.join([mm,mm,mm])] = gulp2mm + gulp3zero
            outdic[''.join([xx,xx,xx])] = gulp2xx + gulp3zero
            
    elif len(atnam) == 3:
        mm, x1, x2 = atnam
        outdic[''.join([mm,x1,x1])] = gulp2 + gulp3mxx
        outdic[''.join([mm,x2,x2])] = gulp2 + gulp3mxx
        outdic[''.join([mm,x1,x2])] = gulp2 + gulp3mx1x2
        outdic[''.join([x1,mm,mm])] = gulp2 + gulp3xmm    
        outdic[''.join([x2,mm,mm])] = gulp2 + gulp3xmm    
        # Add other 2body interactions if Kandemir
        if partype[0] == 'K' or partype[0] == 'k':
            outdic[''.join([mm,mm,mm])] = gulp2mm + gulp3zero
            outdic[''.join([x1,x1,x1])] = gulp2xx + gulp3zero
            outdic[''.join([x2,x2,x2])] = gulp2xx + gulp3zero
    else:
        print('Error in gulpdic: Unexpected number of atoms')
        sys.exit()
        
    return outdic

    
def gulp2lammps(parGulp):
    '''
    Converts GULP SW parameters into the format for LAMMPS.
    
    The transformation is as follows:
        p=4; q=0; epsilon=1;
        A_L = A; sigma = rho; lambda = K
        B_L = B/(rho**4); a = rmax/rho
        gamma = rho12/rho (which is equal to 1 for Jiang's parameterization 
                           but not for Kandemir where rho is different)
        
    INPUT:
        parGulp: Parameter list for SW in GULP
                 (A, B, rho, rmax, K, theta0, rho12, rho13, rmax12, rmax13)
                 
    OUTPUT:
        parLammps: Parameter list for SW in LAMMPS
                  (epsilon, A_L, B_L, sigma, p, q, a, lambda, gamma, costheta0)
    '''
    A, rho, B, rmax, K, theta0, rho12, rho13, rmax12, rmax13 = parGulp
    
    # Transformations
    p=4; q=0; epsilon=1; tol = 0.0
    A_L = A; sigma = rho; lmbda = K
    B_L = B/(rho**4); a = rmax/rho
    gamma = rho12/rho     
    costheta0 = np.cos(theta0*np.pi/180.)
    
    return [epsilon, sigma, a, lmbda, gamma, costheta0, A_L, B_L, p, q, tol]


def genLAMMPSfile(fname, atomName, dicGulp, partype='Jiang2017'):
    '''
    Currently for Jiang type parameterization only
    '''
    # List of possible metal and chalcogen atoms

    
    # Check if input atoms belong to known m and x
    Mlist = []; Xlist = []
    for atom in atomName:
        # Check for name and possible appended name, i.e. "Mo" & "Mo1" both are metals
        if atom in METAL or atom[:-1] in METAL:
            Mlist.append(atom)
            
        elif atom in CHALCOGEN or atom[:-1] in CHALCOGEN:
            Xlist.append(atom)
            
        else:
            print("genLAMMPSfile: Unknown element in atomName")
            print("Define whether %s is metal or chalcogen (M/X)?")
            x = input()
            if x == 'M' or x == 'm':
                Mlist.append(atom)
            elif x == 'X' or x == 'x':
                Xlist.append(atom)
            else:
                print("Wrong input!!")
                sys.exit()
    
    
    # Header for the file
    header = [ "# Generated on %s by Kemal Atalar, using %s parameters"%(date.today(),partype),
             "",
             "# The Stillinger-Weber parameters, for transition-metal dichalcogenide (TMD) lateral heterostructures",
             "# M = " + ", ".join(sorted(set(Mlist))) + "; X = " + ", ".join(sorted(set(Xlist))), 
             "",
             "# these entries are in LAMMPS \"metal\" units:",
             "#   epsilon = eV; sigma = Angstroms",
             "#   other quantities are unitless",
             "",
             "# format of a single entry (one or more lines):",
             "#   element 1, element 2, element 3,",
             "#   epsilon, sigma, a, lambda, gamma, costheta0, A, B, p, q, tol",
             ""
            ]
    
    x1x2_explain = [ "# Numerically labelled chalcogens indicate "
                    ,"# the top and bottom sub-layers within a single layer",
                    ""
            ]
    
    if any([i[-1].isnumeric() for i in atomName]):
        [header.insert(-8,line) for line in x1x2_explain]
        
    ############ SW parameters #############
    
    # Zero list - for interactions that is not described in dicGulp
    zerolist = ["# zero terms"]
    zeroPar = [0] + [1]*7 + [4,0,0.0]
    
    # Initiate all the other parameter groups
    mxx = ["# M-X-X terms"]; 
    xmm = ["# X-M-M terms"]; 
    mx1x2 = ["# M-X1-X2 terms"]
    
    # Iterate over all possible combinations
    comblist = []
    for comb in itertools.combinations(atomName*3, 3):
        comblist.append(comb)
    comblist = set(comblist)
    
    for comb in comblist:
        combstr = ''.join(comb)
        # Check if interaction is defined
        if combstr not in dicGulp:
            zerolist.append(PFORMAT.format(*comb,*zeroPar))
    
        # If it is defined, assign to lists
        else:
            # Convert GULP parameters to LAMMPS format
            lammpsPar = gulp2lammps(dicGulp[combstr])
            line = PFORMAT.format(*comb,*lammpsPar)
            # Check comb type
            if comb[0] in Xlist:
                xmm.append(line)
            else:
                if comb[1] == comb[2]:
                    mxx.append(line)
                else:
                    mx1x2.append(line)
                
    with open(fname, 'w+') as f:
        # Write the header of the potential file
        for line in header + sorted(mxx) + sorted(xmm) + sorted(mx1x2) + sorted(zerolist):
            f.write(line)
            f.write('\n')
    
    return 1


def genMultiple(inpfList,outfname,partype='Jiang2017'):
    '''
    Function to automatically generate  the potential file given the location
    of GULP input files
    '''
    dicGulp={}
    atomName = []
    # Iterate over GULP input files
    for fgulp in inpfList:
        #atName, [gulp2, gulp3_mxx, gulp3_xmm, gulp3_mx1x2] = readGulpInput(fgulp)
        atName, gulpPar = readGulpInput(fgulp)
        atomName += atName
        
        if len(atName) == 2:
            gulp2, gulp3_mxx, gulp3_xmm = gulpPar
            dicGulp.update(gulpDic(atName, gulp2, gulp3_mxx, gulp3_xmm))
        elif len(atName) == 3:
            gulp2, gulp3_mxx, gulp3_xmm, gulp3_mx1x2 = gulpPar
            dicGulp.update(gulpDic(atName, gulp2, gulp3_mxx, gulp3_xmm, gulp3_mx1x2))
    
    genLAMMPSfile(outfname, atomName, dicGulp, partype=partype)
    
    return 1
        
if __name__ == "__main__":
    
    #case = 'MoS2-Jiang'
    #case = 'NbSe2'
    case = 'NbSe2-read' # Automatically read parameters from GULP file
    #case = 'metallic-read'
    
    if case == 'MoS2-Jiang':
        # GULP SW Parameters
        gulp2 = [6.918, 1.252, 17.771, 3.16] # 2-body: A, rho, B, rmax
        # 3-body: K, theta0, rho12, rho13, rmax12, rmax 13
        gulp3_mxx = [67.883, 81.788, 1.252, 1.252, 3.16, 3.16] 
        gulp3_xmm = [62.449, 81.788, 1.252, 1.252, 3.16, 3.16]
        
        at_name = ['Mo','S']
        
        # Old implementation
        #gulpmxx = gulp2 + gulp3_mxx
        #gulpxmm = gulp2 + gulp3_xmm
        #dicGulp = {'MoSS':gulpmxx, 'SMoMo':gulpxmm}
        #lmp_mxx = gulp2lammps(gulpmxx)
        #print(PFORMAT.format('Mo','S','S',*lmp_mxx))
        
        dicGulp = gulpDic(at_name, gulp2, gulp3_mxx, gulp3_xmm)
        
        genLAMMPSfile('MoS2.sw', at_name, dicGulp, partype='Jiang2017')

    
    elif case == 'NbSe2':
        # GULP SW Parameters
        gulp2 = [6.942, 1.138, 22.849, 3.460] # 2-body: A, rho, B, rmax
        # 3-body: K, theta0, rho12, rho13, rmax12, rmax 13
        gulp3_mxx = [34.409, 83.129, 1.138, 1.138, 3.46, 3.46] 
        gulp3_mx1x2 = [34.973, 79.990, 1.138, 1.138, 3.46, 3.46] 
        gulp3_xmm = [34.409, 83.129, 1.138, 1.138, 3.46, 3.46]       
        
        at_name = ['Nb','Se1', 'Se2']
        
        dicGulp = gulpDic(at_name, gulp2, gulp3_mxx, gulp3_xmm, gulp3_mx1x2)
    
        genLAMMPSfile('NbSe2.sw', at_name, dicGulp, partype='Jiang2017')
    
    elif case == 'NbSe2-read':
        
        # Gulp input file
        # These files for the Jiang parameterization can be found in
        # supplement.zip at http://jiangjinwu.org/sw
        fgulp = '/home/ka5118/Downloads/supplement/sw_gulp/h-nbse2.gin'
        
        atName, [gulp2, gulp3_mxx, gulp3_xmm, gulp3_mx1x2] = readGulpInput(fgulp)
 
        dicGulp = gulpDic(atName, gulp2, gulp3_mxx, gulp3_xmm, gulp3_mx1x2)
        
        genLAMMPSfile('NbSe2.sw', atName, dicGulp, partype='Jiang2017')
       
    elif case == 'metallic-read':
        
        # List of input files
        drc = '/home/ka5118/Downloads/supplement/sw_gulp/'
        mat = ['h-nbse2','h-nbs2','h-tas2','h-tase2']
        list_inp = [drc + i + '.gin' for i in mat]
        
        outfnam = 'metallicTMD.sw'
        
        genMultiple(list_inp, outfnam)
       
        
    
    
    
    
    