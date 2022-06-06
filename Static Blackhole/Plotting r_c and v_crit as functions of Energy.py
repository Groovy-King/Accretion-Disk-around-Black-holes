import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import fsolve 

# Defining the relevant constants
r0 = 1
alpha = r0/2
Gamma = 1.50
E = np.linspace(0, 10, 1001)

#Determining the array of critical radii
r_crit = np.array([])
for e in E:
    crit = lambda rc: ((2/(Gamma - 1) + 1/2) * (rc*r0)/(4*(rc - r0)**2) - rc/(2*(rc - r0)) - e)[0]
    r_crit = np.append(r_crit, fsolve(crit, 1.5))

#Determining the speed of sound/gas at sonic point
def Sonic_speed(rc):
    return np.sqrt(rc*r0/(4*(rc - r0)**2))

#Determining the array of critical velocities
v_critical= Sonic_speed(r_crit)

#Plotting r_critical vs E
plt.plot(E, r_crit)
plt.grid()
plt.xlabel("Energy of the system")
plt.ylabel("Critical radius (in r0 units)")
plt.title("Plot of critical radius of the system as a function of Energy")
plt.show()

#Plotting v_critical vs E
plt.plot(E, v_critical)
plt.plot(max(E[v_critical <= 1]), 1, marker='o', markersize=10, markeredgecolor='red', markerfacecolor='red')
plt.grid()
plt.xlabel("Energy of the system")
plt.ylabel("Critical velocity (assuming c = 1)")
plt.title("Plot of velocity at sonic point as a function of Energy")
plt.show()

print(f"The energy at which the velocity at the sonic point equals the speed of light is: {max(E[v_critical <= 1])}")