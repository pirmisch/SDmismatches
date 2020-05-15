# -*- coding: utf-8 -*-
"""
demo of an incumbent mismatch prediction including data

If you use parts of this code please cite:
"""

from SDmismatches import SDMMrates
import numpy as np
import matplotlib.pyplot as plt

#select mode:
mode = "singleIncumbent"

#strand concentrations
c = 10E-9         #10nM

#number of displacement positions
N = 17

#define energy parameters
dGBP = 2.52       #base pairing energy

#define time constant
kbm = 36E3        #branch migration rate constant

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
toe3Samples = np.zeros((sampleSize,N))
toe4Samples = np.zeros((sampleSize,N))
toe5Samples = np.zeros((sampleSize,N))

#initiate produce parameter samples
parameterSamples = np.random.multivariate_normal(popt, pcov,sampleSize)

#calculate MM position dependence for each set of parameters
i=0
for pi in parameterSamples:
    
    for mismatchPosition in range(1,N):
        
        toe3Samples[i,mismatchPosition] = SDMMrates(mismatchPosition,3,N,dGBP,*pi,kbm,c,mode)
        toe4Samples[i,mismatchPosition] = SDMMrates(mismatchPosition,4,N,dGBP,*pi,kbm,c,mode)
        toe5Samples[i,mismatchPosition] = SDMMrates(mismatchPosition,5,N,dGBP,*pi,kbm,c,mode)

    #calculate relative rate
    toe3Samples[i,:] = toe3Samples[i,:] / SDMMrates(0,3,N,dGBP,*pi,kbm,c,mode)
    toe4Samples[i,:] = toe4Samples[i,:] / SDMMrates(0,4,N,dGBP,*pi,kbm,c,mode)
    toe5Samples[i,:] = toe5Samples[i,:] / SDMMrates(0,5,N,dGBP,*pi,kbm,c,mode)    
    
    i = i + 1


#form array of samples
toe3Array = np.asarray(toe3Samples)
toe4Array = np.asarray(toe4Samples)
toe5Array = np.asarray(toe5Samples)


#calculate upper and lower percentiles

toe3lower = np.percentile(toe3Array, 2.5, axis=0)
toe3upper = np.percentile(toe3Array, 97.5, axis=0)

toe4lower = np.percentile(toe4Array, 2.5, axis=0)
toe4upper = np.percentile(toe4Array, 97.5, axis=0)

toe5lower = np.percentile(toe5Array, 2.5, axis=0)
toe5upper = np.percentile(toe5Array, 97.5, axis=0)


##############################################################################
###create optimal solution####################################################

SDrates3 = np.zeros(N)
SDrates4 = np.zeros(N)
SDrates4long = np.zeros(22)
SDrates5 = np.zeros(N)

for mismatchPosition in range(N):
    SDrates3[mismatchPosition] = SDMMrates(mismatchPosition,3,N,dGBP,*popt,kbm,c,mode)
    SDrates4[mismatchPosition] = SDMMrates(mismatchPosition,4,N,dGBP,*popt,kbm,c,mode)
    SDrates5[mismatchPosition] = SDMMrates(mismatchPosition,5,N,dGBP,*popt,kbm,c,mode)
  
for mismatchPosition in range(22):
    SDrates4long[mismatchPosition] = SDMMrates(mismatchPosition,4,22,dGBP,*popt,kbm,c,mode)
##############################################################################
###create plot################################################################
positions=np.arange(1, N, 1)
    
plt.fill_between(positions, toe3lower[1:], toe3upper[1:],alpha=0.5)
plt.fill_between(positions, toe4lower[1:], toe4upper[1:],alpha=0.5)
plt.fill_between(positions, toe5lower[1:], toe5upper[1:],alpha=0.5)

#restart color cycle
plt.gca().set_prop_cycle(None)
 
plt.semilogy(positions,SDrates3[1:]/SDrates3[0])
plt.semilogy(positions,SDrates4[1:]/SDrates4[0])
plt.semilogy(positions,SDrates5[1:]/SDrates5[0])

positions2=np.arange(1, 22, 1)
plt.semilogy(positions2,SDrates4long[1:]/SDrates4long[0],'--')

plt.ylabel("Relative displacement rate")
plt.xlabel("Position of incumbent mismatch")

##############################################################################
###add data###################################################################
#restart color cycle
plt.gca().set_prop_cycle(None)

data3 = np.loadtxt("T3.txt")
data4 = np.loadtxt("T4.txt")
data5 = np.loadtxt("T5.txt")
data4_long = np.loadtxt("T4_long.txt")

plt.errorbar(data3[1:,0],data3[1:,1]/data3[0,1],yerr=data3[1:,1]/data3[0,1]*yerror,fmt='<')
plt.errorbar(data4[1:,0],data4[1:,1]/data4[0,1],yerr=data4[1:,1]/data4[0,1]*yerror,fmt='o')
plt.errorbar(data5[1:,0],data5[1:,1]/data5[0,1],yerr=data5[1:,1]/data5[0,1]*yerror,fmt='D')
plt.errorbar(data4_long[1:,0],data4_long[1:,1]/data4_long[0,1],yerr=data4_long[1:,1]/data4_long[0,1]*yerror,mfc='none',fmt='o')

