# -*- coding: utf-8 -*-
"""
this program integrates a set of differential equations for
a simple reaction network

uses calculated rate constants to predict kinetics

If you use parts of this code please cite:

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from SDmismatches import SDMMrates

#define used functions
def diff(y, t):
    I1,I2,OT,I1OT,I2OT,O,I1T,I2T,I2I1T,I1I2T = y      # unpack current values of y
    derivs = [-k11p*I1*OT  + k11m*I1OT   - k31p*I1*I2T + k31m*I1I2T + k32*I2I1T  ,  #I1
              -k21p*I2*OT  + k21m*I2OT   - k31p*I2*I1T + k31m*I2I1T + k32*I1I2T  ,  #I2
              -k11p*I1*OT  + k11m*I1OT   - k21p*I2*OT  + k21m*I2OT               ,  #OT
              +k11p*I1*OT  - k11m*I1OT   - k12*I1OT                              ,  #I1OT
              +k21p*I2*OT  - k21m*I2OT   - k22*I2OT                              ,  #I2OT
              +k12*I1OT    + k22*I2OT                                            ,  #O
              +k12*I1OT    - k31p*I2*I1T + k31m*I2I1T  + k32*I1I2T               ,  #I1T
              +k22*I2OT    - k31p*I1*I2T + k31m*I1I2T  + k32*I2I1T               ,  #I2T
              +k31p*I2*I1T - k31m*I2I1T  - k32*I2I1T                             ,  #I2I1T
              +k31p*I1*I2T - k31m*I1I2T  - k32*I1I2T                                #I1I2T
            ]
    return derivs

 
#define energy parameters
dGBP    = 2.52       #base pairing energy
dGp     = 3.5        #branch migration initiation energy
dGBM    = 7.4        #branch migration penalty
dGassoc = 2.5        #association energy
dGMM    = 9.5        #mismatch penalty

#define time constant
kbm = 36E3        #branch migration rate constant
    
# Initial concentrations
I1    = 40.0E-9
I2    = 20.0E-9
OT    = 10.0E-9
I1OT  = 0.0
I2OT  = 0.0
O     = 0.0
I1T   = 0.0
I2T   = 0.0
I2I1T = 0.0
I1I2T = 0.0

c=1.0E-9 #can be arbitrary, since k1 and k2 are concentration independent

k11p , k11m , k12 = SDMMrates(18,4,20,dGBP,dGp,dGBM,dGassoc,dGMM,kbm,c,"singleIncumbent",output=3)
k21p , k21m , k22 = SDMMrates(3,4,20,dGBP,dGp,dGBM,dGassoc,dGMM,kbm,c,"singleIncumbent",output=3)
k31p , k31m , k32 = SDMMrates(0,4,24,dGBP,dGp,dGBM,dGassoc,dGMM,kbm,c,"singleIncumbent",output=3,T2=4)


#starting concentration
c0 = [I1,I2,OT,I1OT,I2OT,O,I1T,I2T,I2I1T,I1I2T]


xdata=np.linspace(0,2400,2400)

#solve ODE
y1opt=odeint(diff, c0, xdata, rtol=1e-12, atol=1e-12)

##############################################################################
###add data###################################################################

plt.plot(xdata/60,y1opt[:,5]/1E-8)
plt.plot(xdata/60,y1opt[:,6]/1E-8)
plt.plot(xdata/60,y1opt[:,7]/1E-8)

plt.xlim( -1, 40 )
plt.ylim( -0.1, 1.1 )


plt.xlabel("Time (min)")
plt.ylabel("Normalized product")