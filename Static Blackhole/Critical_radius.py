import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import fsolve

# Defining the relevant constants
r0 = 1
alpha = r0/2
Gamma = 1.50

def Critical_radius(E):
    #Determing the radius of sonic point
    def Crit(rc):
        out = (2/(Gamma - 1) + 1/2) * (rc*r0)/(4*(rc - r0)**2) - rc/(2*(rc - r0)) - E
        return out

    rc = fsolve(Crit, 2)
    print(f"The value of r at the critical point is: {rc[0]} (in r0 units)")
    
    r1 = np.linspace(1.5, 5, 1000)
    plt.plot(r1, Crit(r1))
    plt.grid()
    plt.xlabel("Distance from Center of Black hole")
    plt.ylabel("Value of 2/r - GM/(r-r0)^2 when vr = cs")
    plt.title("Plot of Numerator of RHS at v_r = c_s to determine value of r_c")
    plt.show()
    return rc