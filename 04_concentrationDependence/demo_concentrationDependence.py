# -*- coding: utf-8 -*-
"""
demo how the strand displacement is separated into 
second- and first order contributions

If you use parts of this code please cite:
    
"""

from SDmismatches import SDMMrates
import numpy as np
import matplotlib.pyplot as plt

#select mode:
#mode="singleInvader"
#mode="doubleInvader"
mode="singleIncumbent"

#define energy parameters
dGBP    = 2.52       #base pairing energy
dGp     = 3.5        #branch migration initiation energy
dGBM    = 7.4        #branch migration penalty
dGassoc = 2.5        #association energy
dGMM    = 9.5        #mismatch penalty

#define time constant
kbm = 36E3        #branch migration rate constant

##############################################################################
###single value output########################################################
N                = 17  #number of displacement positions
toeholdLength    = 4
mismatchPosition = 5   #position 0 results in no mismatch

#initiate data vector
dataPoints = 100
keff       = np.zeros(dataPoints)
kfirst     = np.zeros(dataPoints)
ksecond    = np.zeros(dataPoints)

c = np.logspace(-6,-1,dataPoints)

for i in range(dataPoints):
    keff[i] , ksecond[i] , kfirst[i] , cThresh = SDMMrates(mismatchPosition,toeholdLength,N,dGBP,dGp,dGBM,dGassoc,dGMM,kbm,c[i],mode,output=2)
 

plt.loglog(c,keff)
plt.loglog(c,kfirst,'k--')
plt.loglog(c,ksecond*c,'k--')

plt.axvline(x=cThresh,color='grey',linestyle='dotted')

plt.ylim(np.min(keff)*0.9,np.max(keff)*2)

plt.xlabel("Concentration (M)")
plt.ylabel("Displacement rate (1/s)")
