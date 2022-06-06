from ast import For
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import fsolve
from Blackhole_Model import *

# Defining the relevant constants
r0 = 1
Gamma = 1.50

def Critical_radius_Rotating_Blackhole(E, L, a):
    #Determing the radius of sonic point
    def Crit(rc):
        out = (2*Gamma/(Gamma - 1)) * ((L**2/rc**3 - Force(rc, a)) / (dFdr(rc, a) / Force(rc, a) - 3/rc)) + Potential(rc, a) + L**2/(2*rc**2) - E
        return out

    rc = fsolve(Crit, 2.5)
    print(f"The radius of the critical point is: {rc[0]} (in r0 units)")
    return rc