# -*- coding: utf-8 -*-
"""
demo of the function for:
    
**mean reaction times for single incumbent mismatches**

If you use parts of this code please cite:
"""

from SDmismatches import SDMMrates
import numpy as np
import matplotlib.pyplot as plt

#select mode:
#mode = "singleInvader"
#mode = "doubleInvader"
mode = "singleIncumbent"

#strand concentrations
c = 10E-9         #10nM

#number of displacement positions
N = 22

#define energy parameters
dGBP = 2.52       #base pairing energy
dGp = 3.5        #branch migration initiation energy
dGBM = 7.4        #branch migration penalty
dGassoc = 2.5     #association energy
dGMM = 9.5        #mismatch penalty

#define time constant
kbm = 36E3        #branch migration rate constant

##############################################################################
###single value output########################################################
toeholdLength = 4
mismatchPosition = 5 #position 0 results in no mismatch

SDrate = SDMMrates(mismatchPosition,toeholdLength,N,dGBP,dGp,dGBM,dGassoc,dGMM,kbm,c,mode)

print('mean reaction time: {:.1f} s (toehold length: {:.0f} mismatch position: {:.0f})'.format(1/SDrate,toeholdLength,mismatchPosition))


##############################################################################
###create mismatch position dependence #######################################
position = np.arange(1, N, 1)

SDrates = np.zeros(N)


for mismatchPosition in range(N):
    SDrates[mismatchPosition] = SDMMrates(mismatchPosition,toeholdLength,N,dGBP,dGp,dGBM,dGassoc,dGMM,kbm,c,mode)
    
plt.semilogy(position,SDrates[1:N]/SDrates[0])
    
plt.ylabel("Relative displacement rate")
plt.xlabel("Position of incumbent mismatch")