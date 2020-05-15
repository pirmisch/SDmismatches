# -*- coding: utf-8 -*-
"""
demonstrates a fit to single invader-target mismatch data
obtains energy parameters

If you use parts of this code please cite:
    
"""

from SDmismatches import SDMMrates
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def fitFunction(x,param1,param2,param3,param4):
    
    #unpack position vector
    x1 = x[0:size1]
    x2 = x[size1:size2end]
    x3 = x[size2end:size3end]
    x4 = x[size3end:size4end]
    x5 = x[size4end:size5end]
    x6 = x[size5end:size6end]

    #create rate vectors
    y1 = np.zeros(size1)
    y2 = np.zeros(size2)
    y3 = np.zeros(size3)
    y4 = np.zeros(size4)
    y5 = np.zeros(size5)
    y6 = np.zeros(size6)
    
    #get rates
    for index,item in enumerate(x1):
        y1[index] = SDMMrates(int(item),10,N,dGBP,param1,param2,param3,param4,kbm,c,mode)
    for index,item in enumerate(x2):        
        y2[index] = SDMMrates(int(item),9,N,dGBP,param1,param2,param3,param4,kbm,c,mode)
    for index,item in enumerate(x3):   
        y3[index] = SDMMrates(int(item),8,N,dGBP,param1,param2,param3,param4,kbm,c,mode)
    for index,item in enumerate(x4):   
        y4[index] = SDMMrates(int(item),7,N,dGBP,param1,param2,param3,param4,kbm,c,mode)
    for index,item in enumerate(x5):   
        y5[index] = SDMMrates(int(item),6,N,dGBP,param1,param2,param3,param4,kbm,c,mode)
    for index,item in enumerate(x6):   
        y6[index] = SDMMrates(int(item),5,N,dGBP,param1,param2,param3,param4,kbm,c,mode)    
    
    #calculate relative rates
    y1 = y1/y1[0]
    y2 = y2/y2[0]
    y3 = y3/y3[0]
    y4 = y4/y4[0]
    y5 = y5/y5[0]
    y6 = y6/y6[0]

    #recombine rate vector   
    out = np.concatenate((y1,y2,y3,y4,y5,y6))
    return out

#select mode:
mode = "singleInvader"

#strand concentrations
c = 10E-9         #10nM

#number of displacement positions
N = 17

#obtain base pairing energy from DINAmelt, NUPACK,...
dGBP = 2.52      #base pairing energy



#define time constant
kbm = 36E3        
#can be arbitrary in case of a fit to relative displacement rates

##############################################################################
###load data##################################################################

dataT10 = np.loadtxt("T10.txt")
dataT09 = np.loadtxt("T9.txt")
dataT08 = np.loadtxt("T8.txt")
dataT07 = np.loadtxt("T7.txt")
dataT06 = np.loadtxt("T6.txt")
dataT05 = np.loadtxt("T5.txt")


xdataT10 = dataT10[:,0]
ydataT10 = dataT10[:,1] / dataT10[0,1]

xdataT09 = dataT09[:,0]
ydataT09 = dataT09[:,1] / dataT09[0,1]

xdataT08 = dataT08[:,0]
ydataT08 = dataT08[:,1] / dataT08[0,1]

xdataT07 =dataT07[:,0]
ydataT07 =dataT07[:,1] / dataT07[0,1]

xdataT06 = dataT06[:,0]
ydataT06 = dataT06[:,1] / dataT06[0,1]

xdataT05 = dataT05[:,0]
ydataT05 = dataT05[:,1] / dataT05[0,1]

#get sizes of data vectors
size1 = np.size(xdataT10)
size2 = np.size(xdataT09)
size3 = np.size(xdataT08)
size4 = np.size(xdataT07)
size5 = np.size(xdataT06)
size6 = np.size(xdataT05)

size2end = size1    + size2
size3end = size2end + size3
size4end = size3end + size4
size5end = size4end + size5
size6end = size5end + size6

#create a combined data vectors
xdata = np.concatenate((xdataT10, xdataT09, xdataT08, xdataT07, xdataT06, xdataT05), axis=0)
ydata = np.concatenate((ydataT10, ydataT09, ydataT08, ydataT07, ydataT06, ydataT05), axis=0)

#define errors
error  = 0.35
yerror = error * ydata 

##############################################################################
###perform fit################################################################

#define start parameters for fit
dGp     = 3           #branch migration initiation energy
dGBM    = 8          #branch migration penalty
dGassoc = 3       #association energy
dGMM    = 9       #mismatch penalty

pstart = (dGp,dGBM,dGassoc,dGMM)

#perform weighted fit
popt,pcov = curve_fit(fitFunction, xdata, ydata, sigma = yerror, absolute_sigma=True, p0=pstart)

#print(popt) #used for prediction
#print(pcov) #used for prediction

#calculate errors
pErr = np.sqrt(np.diag(pcov))


print('dGp ={:.1f}, error: {:.1f}'.format(popt[0],pErr[0]))
print('dGBM ={:.1f}, error: {:.1f}'.format(popt[1],pErr[1]))
print('dGassoc ={:.1f}, error: {:.1f}'.format(popt[2],pErr[2]))
print('dGMM ={:.1f}, error: {:.1f}'.format(popt[3],pErr[3]))

##############################################################################
###create mismatch position dependence (Figure 3A)############################

position = np.arange(1, N, 1)


for toeholdLength in range(5,11):
    SDrates = np.zeros(N)
    

    for mismatchPosition in range(N):
        SDrates[mismatchPosition] = SDMMrates(mismatchPosition,toeholdLength,N,dGBP,*popt,kbm,c,mode)
        
    plt.semilogy(position,SDrates[1:N]/SDrates[0])
    


##############################################################################
###plot measurements (Figure 3A)############################
    
#restart color cycle
plt.gca().set_prop_cycle(None)

#plot data    
plt.errorbar(xdataT05[1:],ydataT05[1:],yerr=ydataT05[1:]*error,fmt='v')
plt.errorbar(xdataT06[1:],ydataT06[1:],yerr=ydataT06[1:]*error,fmt='D') 
plt.errorbar(xdataT07[1:],ydataT07[1:],yerr=ydataT07[1:]*error,fmt='>')
plt.errorbar(xdataT08[1:],ydataT08[1:],yerr=ydataT08[1:]*error,fmt='o')
plt.errorbar(xdataT09[1:],ydataT09[1:],yerr=ydataT09[1:]*error,fmt='^')
plt.errorbar(xdataT10[1:],ydataT10[1:],yerr=ydataT10[1:]*error,fmt='s')





    
plt.ylabel("Relative displacement rate")
plt.xlabel("Position of mismatch")

