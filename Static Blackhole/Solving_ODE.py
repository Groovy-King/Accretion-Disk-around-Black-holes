import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import fsolve
from scipy.integrate import odeint
from sympy import rsolve

# Defining the relevant constants
r0 = 1
alpha = r0/2
Gamma = 1.50

def Solve_ODE(rc, vr_crit, ini_guess):
    cs_crit = vr_crit
    def der(S):
        dcdr, dvdr = S
        num1 = ((Gamma - 1)/2) * (4*vr_crit*dvdr/rc - 2*vr_crit**2/rc**2 + r0/(rc - r0)**3)
        den1 = dcdr - 2*dvdr + dcdr
        out1 = dcdr - num1/den1
        num2 = (4*cs_crit*dcdr/rc - 2*cs_crit**2/rc**2 + r0/(rc - r0)**3)
        den2 = dvdr - 2*dcdr + dvdr
        out2 = dvdr - num2/den2
        return [out1[0], out2[0]]
    
    der_crit = fsolve(der, ini_guess)
    
    ## Returns the derivative dCs/dr
    def dCdr_accretion(r, vr, cs):
        if (vr != cs):
            return (Gamma - 1)/2 * (2 * vr**2 / r - alpha/(r - r0)**2) / (cs - vr**2/cs)
        else:
            derivative = der_crit[0]
            return derivative

    ## Returns the derivative dVr/dr
    def dVdr_accretion(r, vr, cs):
        if (vr != cs):
            return (2 * cs**2 / r - alpha/(r - r0)**2) / (vr - cs**2/vr)
        else:
            derivative = der_crit[1]
            return derivative
        
    ## Returns the array of Derivatives
    def dSdr_accretion(r, S):
        cs, vr = S 
        return np.array([dCdr_accretion(r, vr, cs), dVdr_accretion(r, vr, cs)])

    #Solving the ODE above rc
    r_solve_above = np.arange(rc, rc + 5, 0.001)
    solution_above = odeint(dSdr_accretion, y0 = np.array([cs_crit[0], vr_crit[0]]), t = r_solve_above, tfirst = True)

    #Solving the ODE Below rc
    r_solve_below = np.arange(rc, r0 + 0.01, -0.001)
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
    if der_crit[1] < 0:
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
    print(f"The radius of breakdown for " + str + " is: {min(r_solve[x <= 1])}")       
    
    return final_solution, r_solve