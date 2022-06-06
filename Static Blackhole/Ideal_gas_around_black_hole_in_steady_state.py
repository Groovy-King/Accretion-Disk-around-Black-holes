import numpy as np
import matplotlib.pyplot as plt

from Critical_radius import Critical_radius
from Solving_ODE import Solve_ODE


# Defining the relevant constants
r0 = 1
alpha = r0/2
Gamma = 1.50
E = float(input("Please input the energy to be used in the calculation:"))

rc = Critical_radius(E)

#Determining the speed of sound/gas at sonic point
def Sonic_speed(rc):
    return np.sqrt(rc*r0/(4*(rc - r0)**2))

cs_crit = Sonic_speed(rc)
vr_crit = cs_crit
print(f"The velocity of gas at the sonic point is {vr_crit[0]}")

#Solving the ODE system for Accretion
ini_guess = (-0.1, -0.2)
sol_accretion, r_solve = Solve_ODE(rc, vr_crit, ini_guess)


## We repeat the same process for Wind
ini_guess = (-0.2, 0.2)
sol_wind, r_solve = Solve_ODE(rc, vr_crit, ini_guess)


plt.plot(r_solve, sol_accretion, label="Accretion")
plt.plot(r_solve, sol_wind, label="Wind")
plt.plot(rc, 1, marker='o', markersize=10, markeredgecolor='red', markerfacecolor='green')
plt.grid()
plt.xlabel("The distance from center of black hole (in r0 units)")
plt.ylabel("Ratio of speed of gas to speed of sound")
plt.title(f"Plot for v_r/c_s with Energy, E = {E}, with r_critical = {round(rc[0], 2)} r0")
plt.legend()
plt.show()