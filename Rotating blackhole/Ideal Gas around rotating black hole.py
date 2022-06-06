import numpy as np
import matplotlib.pyplot as plt

from Critical_radius import Critical_radius_Rotating_Blackhole
import Blackhole_Model as bhm
from Solving_ODE import Solve_ODE

# Defining the relevant constants
r0 = 1
Gamma = 1.50
p = np.sqrt(2/(Gamma + 1))
E = float(input("Please input the energy to be used in the calculation: "))
L = float(input("Please input the Angular momentum of the particle: "))
A = float(input("Please input the Angular momentum of the black hole (0 <= A <= 1): "))

rc = Critical_radius_Rotating_Blackhole(E, L, A)

#Determining the velociy at critical point
cs_crit = np.sqrt((Gamma + 1) * (L**2/rc**3 - bhm.Force(rc, A)) / (bhm.dFdr(rc, A)/bhm.Force(rc, A) - 3/rc))
vr_crit = np.sqrt(2 / (Gamma + 1)) * cs_crit

print(f"The velocity of the gas at the critical point is: {vr_crit[0]}")

#Solving the ODE system for Accretion
sol_accretion, r_solve = Solve_ODE(rc, cs_crit, L, A, 1)


## We repeat the same process for Wind
sol_wind, r_solve = Solve_ODE(rc, cs_crit, L, A, -1)


plt.plot(r_solve, sol_accretion, label="Accretion")
plt.plot(r_solve, sol_wind, label="Wind")
plt.plot(rc, p, marker='o', markersize=10, markeredgecolor='red', markerfacecolor='green')
plt.grid()
plt.xlabel("The distance from center of black hole (in r0 units)")
plt.ylabel("Ratio of speed of gas to speed of sound")
plt.title(f"Plot for v_r/c_s with E = {E}, L = {L}, a = {A}, with r_critical = {round(rc[0], 2)} r0")
plt.legend()
plt.show()