import numpy as np
import sympy as sp
from scipy.integrate import quad

#Defining the variables
r, a = sp.symbols('r, a')
num = (r**2 - 2*a*sp.sqrt(r) + a**2)**2
den = r**3 * (sp.sqrt(r)*(r - 2) + a)**2
force_sym = num/den

Force = sp.lambdify((r, a), force_sym)

dFdr = sp.lambdify((r, a), sp.diff(force_sym, r))

dFdr2 = sp.lambdify((r, a), sp.diff(force_sym, r, 2))

#Creating the Potential function
def Potential(r, a):
    return quad(Force, r, np.inf, args=(a))[0]
