# -*- coding: utf-8 -*-
"""
demo of the function for:
    
**mean reaction times for two invader mismatches**

If you use parts of this code please cite:
"""

from SDmismatches import SDMMrates
import numpy as np
import matplotlib.pyplot as plt

#select mode:
#mode = "singleInvader"
mode = "doubleInvader"
#mode = "singleIncumbent"

#strand concentrations
c = 10E-9         #10nM

#number of displacement positions
N = 27

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
toeholdLength = 10
mismatchPosition = 2
mismatchPosition2 = 20

SDrate = SDMMrates(mismatchPosition2,toeholdLength,N,dGBP,dGp,dGBM,dGassoc,dGMM,kbm,c,mode,mismatchPosition)

print('mean reaction time: {:.1f} s (toehold length: {:.0f} mismatch positions: {:.0f} and {:.0f})'.format(1/SDrate,toeholdLength,mismatchPosition,mismatchPosition2))


##############################################################################
###create mismatch position dependence #######################################
position = np.arange(1, N, 1)

SDrates = np.zeros(N)


for mismatchPosition2 in range(N):
    SDrates[mismatchPosition2] = SDMMrates(mismatchPosition2,toeholdLength,N,dGBP,dGp,dGBM,dGassoc,dGMM,kbm,c,mode,mismatchPosition)
    
perfectMatchingRate = SDMMrates(0,toeholdLength,N,dGBP,dGp,dGBM,dGassoc,dGMM,kbm,c,mode,0)

plt.semilogy(position,SDrates[1:N]/perfectMatchingRate)
    
plt.ylabel("Relative displacement rate")
plt.xlabel("Position of second mismatch")