from cProfile import label
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import sympy as smp

from scipy.optimize import fsolve, root


#Defining Constants
r0 = 1
alpha = 1
Gamma = 1.5
E = float(input("Please enter the energy to be used in the calculation:"))

#Determing the radius of sonic point
def Crit(rc):
    out = (2/(Gamma - 1) + 1/2) * (rc*r0)/(4*(rc - r0)**2) - rc/(2*(rc - r0)) - E
    return out

rc = fsolve(Crit, 2)
print(rc)

#Determining the derivatives at the sonic point
cs_crit = np.sqrt(rc*r0/(4*(rc - r0)**2))
vr_crit = cs_crit
print(vr_crit)

def der(S):
    dcdr, dvdr = S
    num1 = ((Gamma - 1)/2) * (4*vr_crit*dvdr/rc - 2*vr_crit**2/rc**2 + r0/(rc - r0)**3)
    den1 = dcdr - 2*dvdr + dcdr
    out1 = dcdr - num1/den1
    num2 = (4*cs_crit*dcdr/rc - 2*cs_crit**2/rc**2 + r0/(rc - r0)**3)
    den2 = dvdr - 2*dcdr + dvdr
    out2 = dvdr - num2/den2
    return [out1[0], out2[0]]

ini_guess = (-0.1, -0.2)
der_crit = fsolve(der, ini_guess)
print(der_crit)

def test(x):
    return x