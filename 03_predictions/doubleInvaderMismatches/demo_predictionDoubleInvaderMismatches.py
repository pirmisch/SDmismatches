# -*- coding: utf-8 -*-
"""
demo of a double invader mismatch prediction including data

If you use parts of this code please cite:
"""

from SDmismatches import SDMMrates
import numpy as np
import matplotlib.pyplot as plt

#select mode:
mode = "doubleInvader"

#strand concentrations
c = 10E-9         #10nM

#number of displacement positions
N = 27

#define energy parameters
dGBP = 2.52       #base pairing energy

#define time constant
kbm = 36E3        #branch migration rate constant

#define toehold length
toeholdLength = 10

#define yerror
yerror = 0.35

# use optimal parameters and covariance matrix obtained by fitting procedure
popt=[3.54400181 , 7.38636903 , 2.50249644 , 9.50522577]

pcov=[[ 0.05483743 , -0.01913771 ,  0.01330909 , -0.02427002], 
      [-0.01913771 ,  0.04258403 ,  0.00961577 , -0.01345973], 
      [ 0.01330909 ,  0.00961577 ,  0.03751586 , -0.00560454], 
      [-0.02427002 , -0.01345973 , -0.00560454 ,  0.03164691]]

##############################################################################
###create confidence intervals################################################

#define sample size
sampleSize = 1000

#initiate samples
singleMMSamples = np.zeros((sampleSize,N))
doubleMMSamples = np.zeros((sampleSize,N))

#initiate produce parameter samples
parameterSamples = np.random.multivariate_normal(popt, pcov,sampleSize)

#calculate MM position dependence for each set of parameters
i=0
for pi in parameterSamples:
    
    for mismatchPosition in range(1,N):
        
        singleMMSamples[i,mismatchPosition] = SDMMrates(mismatchPosition,toeholdLength,N,dGBP,*pi,kbm,c,mode,0)
        doubleMMSamples[i,mismatchPosition] = SDMMrates(mismatchPosition,toeholdLength,N,dGBP,*pi,kbm,c,mode,2)

    #calculate relative rate
    singleMMSamples[i,:] = singleMMSamples[i,:] / SDMMrates(0,toeholdLength,N,dGBP,*pi,kbm,c,mode,0)
    doubleMMSamples[i,:] = doubleMMSamples[i,:] / SDMMrates(0,toeholdLength,N,dGBP,*pi,kbm,c,mode,0)  
    
    i = i + 1


#form array of samples
singleMMArray = np.asarray(singleMMSamples)
doubleMMArray = np.asarray(doubleMMSamples)

#calculate upper and lower percentiles

singleMMlower = np.percentile(singleMMArray, 2.5, axis=0)
singleMMupper = np.percentile(singleMMArray, 97.5, axis=0)

doubleMMlower = np.percentile(doubleMMArray, 2.5, axis=0)
doubleMMupper = np.percentile(doubleMMArray, 97.5, axis=0)


##############################################################################
###create optimal solution####################################################

SDratessingleMM = np.zeros(N)
SDratesdoubleMM = np.zeros(N)

for mismatchPosition in range(N):
    SDratessingleMM[mismatchPosition] = SDMMrates(mismatchPosition,toeholdLength,N,dGBP,*popt,kbm,c,mode,0)
    SDratesdoubleMM[mismatchPosition] = SDMMrates(mismatchPosition,toeholdLength,N,dGBP,*popt,kbm,c,mode,2)

    
##############################################################################
###create plot################################################################
positions=np.arange(1, N, 1)
    
plt.fill_between(positions, singleMMlower[1:], singleMMupper[1:],alpha=0.5)
plt.fill_between(positions, doubleMMlower[1:], doubleMMupper[1:],alpha=0.5)

#restart color cycle
plt.gca().set_prop_cycle(None)
 
plt.semilogy(positions,SDratessingleMM[1:N]/SDratessingleMM[0])
plt.semilogy(positions,SDratesdoubleMM[1:N]/SDratessingleMM[0])

plt.ylabel("Relative displacement rate")
plt.xlabel("Position of incumbent mismatch")


##############################################################################
###add data###################################################################
#restart color cycle
plt.gca().set_prop_cycle(None)

dataSingle = np.loadtxt("T10_single.txt")
dataDouble = np.loadtxt("T10_double.txt")

plt.errorbar(dataSingle[1:,0],dataSingle[1:,1]/dataSingle[0,1],yerr=dataSingle[1:,1]/dataSingle[0,1]*yerror,fmt='s')
plt.errorbar(dataDouble[:,0],dataDouble[:,1]/dataSingle[0,1],yerr=dataDouble[:,1]/dataSingle[0,1]*yerror,fmt='^')