# -*- coding: utf-8 -*-
"""
function to solve first passage model for:
    -single invader mismatches
    -double invader mismatches
    -single incumbent mismatches

input:
    -mismatch position
    -toehold length
    -displacement length
    -energy parameters dGBP,dGp,dGBM,dGassoc,dGMM
    -branch migration rate
    -concentration
    -mismatch mode (single/double invader, single incumbent)
    -optional second mismatch
    -optional second toehold (handle with care)
    
returns:
    
    if output = 1 (preset)
    -inverse of mean first passage time
    
    if output = 2:
    -inverse of mean first passage time
    -first order rate constant
    -second order rate constant
    -threshold concentration
    
    if output = 3:
    -k1 forward  (second order rate constant for toehold binding)
    -k1 backward (first order rate constant for toehold unbinding)
    -k2          (first order rate constant for full displacement)

If you use parts of this code please cite:

    
"""

import numpy as np

def SDMMrates(MMpos,T,N,dGBP,dGp,dGBM,dGassoc,dGMM,kbm,c,mode,MMpos2=None,T2=0,output=1):
 
#check for input errors
    
    if mode!="singleInvader" and mode!="doubleInvader" and mode!="singleIncumbent":
        print("please select a valid mismatch mode!")
        print("valid modes are:")    
        print("singleInvader,doubleInvader,singleIncumbent")    
        return
    
    if MMpos<0:
        print("Negativ mismatch position!")
        print("please use a position in [0,N)!")
        return
    
    if MMpos>=N:
        print("Mismatch position bigger N!")
        print("please use a position in [0,N)!")
        return
    
#create forward and backward rates for all positions
    
    #initialization
    rates = np.vstack([np.linspace(kbm,kbm,T+N+1),np.linspace(kbm,kbm,T+N+1)])
    
    #calculate base pairing rate
    kbp = kbm * np.exp(dGBM)
    
    #introduce boundaries
    rates[0,0]   = 0    
    rates[1,N+T] = 0

    #initial binding
    kbind = kbp * np.exp(-dGassoc) * c
    
    rates[1,0] = kbind
    rates[0,1] = kbp * np.exp(-dGBP)
            
    #toehold zipping
    for z in range (1,T):
        rates[1,z] = kbp
        rates[0,z+1] = kbp * np.exp(-dGBP)
        
    #start branch migration
    rates[1,T] =rates[1,T] * np.exp(-dGp)
    
#introduce mismatch
    
    if mode=="singleInvader":
        if MMpos>0:
            rates[0,T+MMpos] = kbm * np.exp(dGMM)
            
    if mode=="doubleInvader":
        
        #sort out single mismatch cases
        if MMpos>0 and MMpos2==0:
            rates[0,T+MMpos] = kbm * np.exp(dGMM)               
        if MMpos==0 and MMpos2>0:
            rates[0,T+MMpos2] = kbm * np.exp(dGMM)
        if MMpos==MMpos2:              
            rates[0,T+MMpos2] = kbm * np.exp(dGMM)
            
        #calculate the combined macrostate
        if MMpos>0 and MMpos2>0 and MMpos!=MMpos2:
            Z=np.exp(-2*dGMM) + np.exp(-(dGMM+(dGBP+0.8)*np.abs(MMpos-MMpos2)))
            dGtot=-np.log(Z)
            rates[0,T+MMpos] = kbm * np.exp(dGtot/2)
            rates[0,T+MMpos2] = kbm * np.exp(dGtot/2)

    if mode=="singleIncumbent":
        if MMpos>0:
            
            #in case of incumbent-target mismatch 
            rates[1,T] = rates[1,T] * np.exp(dGp)
            
            #initiate a position dependent energy profile
            pos = np.arange(1, MMpos, 1)
            EnergyProfile = np.zeros(MMpos+1)
            
            #free energy for states after mismatch elimination
            EnergyProfile[MMpos] = dGp + np.log(np.exp(-MMpos*dGBP-dGp) + np.exp(-dGMM))
            
            #free energy for states before mismatch elimination
            EnergyProfile[pos] = EnergyProfile[MMpos] - np.log(np.exp(-(MMpos-pos)*dGBP) + np.exp(-dGMM))

            #get energy difference for all states
            ED=np.diff(EnergyProfile)
            
            #tune all rates according to energy difference
            for i in range(0,MMpos):
                rates[1,T+i] = rates[1,T+i] * np.exp(-ED[i])
            
    #for a second toehold introduce T2 mismatches at the end
    for i in range(T2):
        rates[0,T+N-i] = kbm * np.exp(dGMM)
           
#calculate rate for spontaneous dissociation of incumbent
    rateoff =np.zeros(T+N+1)
    for i in range (1,N+1):
        rateoff[T+N-i] = kbp * np.exp(-dGBP*i)  
    
#Theory first passage time
    
    #initialization of the occupation probability and the flux
    poccup = np.zeros(1+T+N)
    rflux = np.zeros(1+T+N) 
    
    #introduce boundaries    
    poccup[T+N] = 0
    poccup[T+N-1] = 1 / kbm
    rflux[T+N-1] = 1 

        
    #Calculate occupancies from rates
    for n in range(2,T+N+1):
        rflux[T+N-n] = rflux[T+N-n+1] + rateoff[T+N-n+1] * poccup[T+N-n+1]
        poccup[T+N-n] = rates[0,T+N-n+1] / rates[1,T+N-n] * poccup[T+N-n+1] + rflux[T+N-n] / rates[1,T+N-n]
    
    #Calculate mean first passage time
    firstPassageTime = sum(poccup) / rflux[0]

    if output==1:
        return 1 / firstPassageTime
    
    #get first order rate constant
    firstOrderConstant = rflux[0] / sum(poccup[1:]) 
    
    #get second order rate constant
    #positive flux
    jp = rflux[1] / poccup[1]
    
    #negative flux
    jm = kbp * np.exp(-dGBP)
    
    #passge probability
    ppassage = jp / (jm + jp)    
    
    #get second order rate constant
    secondOrderConstant = kbind * ppassage / c
    
    #threshold concentration
    cThresh = firstOrderConstant / secondOrderConstant
  
    if output==2:
        return 1 / firstPassageTime,secondOrderConstant,firstOrderConstant,cThresh
    
    if output==3:
        
        #calculate constants for 3 rate approximation
        k1p = kbind / c * (1-np.exp(-dGBP)) / (1-np.exp(-T*dGBP))
        k2  = firstOrderConstant
        k1m = k1p * k2 / secondOrderConstant - k2
        
        return k1p , k1m , k2
    else: 
        print('Wrong output option, chose 1, 2 or 3!')
        return