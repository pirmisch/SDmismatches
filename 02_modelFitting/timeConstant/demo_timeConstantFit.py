# -*- coding: utf-8 -*-
"""
demonstrates a fit to single invader-target mismatch data
obtains time constant

If you use parts of this code please cite:
    
"""

from SDmismatches import SDMMrates
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def fitFunction(x,param):
    
    #unpack position vector
    x1 = x[0:size1]
    x2 = x[size1:size2end]
    x3 = x[size2end:size3end]

    #create rate vectors
    y1 = np.zeros(size1)
    y2 = np.zeros(size2)
    y3 = np.zeros(size3)
    
    #get rates
    for index,item in enumerate(x1):
        y1[index] = SDMMrates(0,int(item),N,dGBP,dGp,dGBM,dGassoc,dGMM,param,c,mode)
    for index,item in enumerate(x2):        
        y2[index] = SDMMrates(2,int(item),N,dGBP,dGp,dGBM,dGassoc,dGMM,param,c,mode)
    for index,item in enumerate(x3):   
        y3[index] = SDMMrates(11,int(item),N,dGBP,dGp,dGBM,dGassoc,dGMM,param,c,mode) 
    
    #recombine rate vector   
    out = np.concatenate((y1,y2,y3))
    
    return out

#select mode:
mode = "singleInvader"

#strand concentrations
c = 10E-9         #10nM

#number of displacement positions
N = 17

#obtain base pairing energy from DINAmelt, NUPACK,...
dGBP = 2.52       #base pairing energy

#fitted energy parameters
dGp     = 3.5           #branch migration initiation energy
dGBM    = 7.4          #branch migration penalty
dGassoc = 2.5       #association energy
dGMM    = 9.5       #mismatch penalty



      
#can be arbitrary in case of a fit to relative displacement rates

##############################################################################
###load data##################################################################

dataP    = np.loadtxt("perfect.txt")
dataMM2  = np.loadtxt("mm2.txt")
dataMM11 = np.loadtxt("mm11.txt")


xdataP = dataP[:,0]
ydataP = dataP[:,1]*c

xdataMM2 = dataMM2[:,0]
ydataMM2 = dataMM2[:,1]*c

xdataMM11 = dataMM11[:,0]
ydataMM11 = dataMM11[:,1]*c


#get sizes of data vectors
size1 = np.size(xdataP)
size2 = np.size(xdataMM2)
size3 = np.size(xdataMM11)

size2end = size1    + size2
size3end = size2end + size3

#create a combined data vectors
xdata = np.concatenate((xdataP,xdataMM2,xdataMM11), axis=0)
ydata = np.concatenate((ydataP,ydataMM2,ydataMM11), axis=0)

#define errors
error = 0.25
yerror = error*ydata 

##############################################################################
###perform fit################################################################

#define start parameters for fit
kbm = 1E5

pstart = kbm

#perform fit
popt,pcov = curve_fit(fitFunction, xdata, ydata,  p0=pstart)

#calculate errors
pErr = np.sqrt(np.diag(pcov))


print('kbm ={:.1e}, error: {:.0e}'.format(popt[0],pErr[0]))

##############################################################################
###create mismatch position dependence (Figure 3B)############################


SDratesP    = np.zeros(size1)
SDratesMM2  = np.zeros(size2)
SDratesMM11 = np.zeros(size3)



for index,item in enumerate(xdataP):
    SDratesP[index] = SDMMrates(0,int(item),N,dGBP,dGp,dGBM,dGassoc,dGMM,*popt,c,mode)

for index,item in enumerate(xdataMM2):
    SDratesMM2[index] = SDMMrates(2,int(item),N,dGBP,dGp,dGBM,dGassoc,dGMM,*popt,c,mode)

for index,item in enumerate(xdataMM11):
    SDratesMM11[index] = SDMMrates(11,int(item),N,dGBP,dGp,dGBM,dGassoc,dGMM,*popt,c,mode)

    
plt.semilogy(xdataP,SDratesP)
plt.semilogy(xdataMM2,SDratesMM2)
plt.semilogy(xdataMM11,SDratesMM11)


##############################################################################
###plot measurements (Figure 3B)############################    

#restart color cycle
plt.gca().set_prop_cycle(None)

plt.errorbar(xdataP,ydataP,yerr=ydataP*error,fmt='v')
plt.errorbar(xdataMM2,ydataMM2,yerr=ydataMM2*error,fmt='^')
plt.errorbar(xdataMM11,ydataMM11,yerr=ydataMM11*error,fmt='s')
    
plt.ylabel("Displacement rate")
plt.xlabel("Position of mismatch")

