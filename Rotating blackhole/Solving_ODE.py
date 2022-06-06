import numpy as np
import matplotlib.pyplot as plt
import sympy as smp

from scipy.optimize import fsolve
from scipy.integrate import odeint
from sympy import rsolve
from Blackhole_Model import *

# Defining the relevant constants
r0 = 1
Gamma = 1.50
p = np.sqrt(2/(Gamma + 1))

def Solve_ODE(rc, cs_crit, L, A, sgn):
    vr_crit = p * cs_crit
    
    #Solving for the derivatives at the sonic point
    l = 1 + 2 / ((Gamma + 1) * p**2) + 4 * (Gamma - 1)**2 / (p**2 * (Gamma + 1)**2)
    m = 4 * (Gamma - 1) * (3/rc - dFdr(rc, A)/Force(rc, A)) * cs_crit / (p * (Gamma + 1)**2)
    n = cs_crit**2 * ((dFdr(rc, A) / Force(rc, A))**2 - dFdr2(rc, A)/Force(rc, A) - 3/rc**2) / (Gamma + 1) - 3*L**2/rc**4 - dFdr(rc, A) - cs_crit**2 * (Gamma - 1) * (3/rc - dFdr(rc, A)/Force(rc, A))**2 / (Gamma + 1)**2
    dvdr_crit = - sgn * (m + np.sqrt(m**2 - 4*l*n)) / (2*l)
    dcdr_crit = - sgn * (Gamma - 1) * dvdr_crit / ((Gamma + 1) * p)
    
    ## Returns the derivative dCs/dr
    def dCdr(r, vr, cs):
        if (vr != p * cs):
            return (Gamma - 1)/2 * (L**2/r**3 - Force(r, A) + vr**2 * (3/r - dFdr(r, A)/Force(r, A))/2 ) / (cs - (Gamma + 1)/2 * vr**2/cs)
        else:
            derivative = dcdr_crit[0]
            return derivative

    ## Returns the derivative dVr/dr
    def dVdr(r, vr, cs):
        if (vr != p * cs):
            return (L**2/r**3 - Force(r, A) + cs**2 * (3/r - dFdr(r, A)/Force(r, A))/(Gamma + 1) ) / (vr - 2/(Gamma + 1) * cs**2/vr)
        else:
            derivative = dvdr_crit[0]
            return derivative
        
    ## Returns the array of Derivatives
    def dSdr_accretion(r, S):
        cs, vr = S 
        return np.array([dCdr(r, vr, cs), dVdr(r, vr, cs)])

    #Solving the ODE above rc
    r_solve_above = np.arange(rc, rc + 5, 0.001)
    solution_above = odeint(dSdr_accretion, y0 = np.array([cs_crit[0], vr_crit[0]]), t = r_solve_above, tfirst = True)

    #Solving the ODE Below rc
    r_solve_below = np.arange(rc, r0 + 0.05, -0.001)
    solution_below = odeint(dSdr_accretion, y0 = np.array([cs_crit[0], vr_crit[0]]), t = r_solve_below, tfirst = True)

    #Reversing the np arrays for the solution below rc
    r_solve_below = np.flipud(r_solve_below)
    solution_below = np.flipud(solution_below)

    #Combining the arrays for the solution
    r_solve = np.concatenate((r_solve_below, r_solve_above))
    solution = np.concatenate((solution_below, solution_above))

    #Determining v_r and c_s
    v_r = np.array([])
    c_s = np.array([])
    for i in range(len(r_solve)):
        v_r = np.concatenate((v_r, np.array([solution[i][1]])))
        c_s = np.concatenate((c_s, np.array([solution[i][0]])))
        
    #For Accretion
    if dvdr_crit < 0:
        str = "Accretion"
        x = v_r
        
    else:
        str = "Wind"
        x = c_s
    
    #Plotting the final solution
    plt.plot(r_solve, c_s, label = "Speed of sound")
    plt.plot(r_solve, v_r, label = "Speed of gas")
    plt.plot(min(r_solve[x <= 1]), 1, label = "Point when c_s = speed of light", marker='o', markersize=10, markeredgecolor='red', markerfacecolor='blue')
    plt.grid()
    plt.xlabel("The distance from center of black hole (in r0 units)")
    plt.ylabel("Speed of sound/gas (assuming speed of light, c = 1)")
    plt.title("Plot of v_r and c_s vs r for " + str)
    plt.legend()
    plt.show()

    #Plotting the ratio of v_r and c_s
    plt.plot(r_solve, v_r/c_s)
    plt.grid()
    plt.xlabel("The distance from center of black hole (in r0 units)")
    plt.ylabel("Ratio of speed of gas to speed of sound")
    plt.title("Plot of v_r/c_s as a function of r in case of " + str)
    plt.show()

    final_solution = v_r/c_s
    print(f"The radius of breakdown for " + str + f" is: {min(r_solve[x <= 1])}")       
    
    return final_solution, r_solve