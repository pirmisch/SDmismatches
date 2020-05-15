# Modelling DNA-strand displacement reactions in the presence of base-pair mismatches
A simple thermodynamic model, which is capable of quantitatively describing the kinetics of strand displacement reactions in the presence of base-pair mismatches, using a minimal set of parameters.

## SDmismatches.py 

solves first passage model for:
* single invader mismatches
* double invader mismatches
* single incumbent mismatches

### input:
* mismatch position
* toehold length
* displacement length
* energy parameters dGBP,dGp,dGBM,dGassoc,dGMM
* branch migration rate
* concentration
* mismatch mode (single/double invader, single incumbent)
* optional second mismatch
* optional second toehold (handle with care)
    
### returns:

if output = 1 (preset)
*inverse of mean first passage time
    
if output = 2:
*inverse of mean first passage time
*first order rate constant
*second order rate constant
*threshold concentration
    
if output = 3:
*k1 forward  (second order rate constant for toehold binding)
*k1 backward (first order rate constant for toehold unbinding)
*k2          (first order rate constant for full displacement)

## Demo
### 01 basic functions
calculates mean-first-passage-time and plots mismatch position dependence for:
* single invader-target mismatches
* double invader-target mismatches
* single incumbent-target mismatches

### 02 model fitting
uses measured data to perform a fit for:
* energy parameters
* branch migration rate constant

### 03 predictions
compares measured data to a prediction using the parameters obtained in 02 for: 
* incumbent-target mismatches
* double invader-target mismatches

### 04 Concentration dependence (demo for output option 2)
The overall rate constant is seprated into a second- and a first-order contribution:
* concentration dependence


### 05 DNA network (demo for output option 3)
The overall rate constant is seprated into three rate constants to set up a system of ordinary differential equations:
* pulse generation using a DNA network

## References
'Modelling DNA-strand displacement reactions in the presence of base-pair mismatches' Patrick Irmisch, Thomas E. Ouldridge, Ralf Seidel
-link to publication
