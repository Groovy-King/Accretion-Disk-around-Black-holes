import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import sympy as smp

from scipy.optimize import fsolve
from scipy.integrate import quad, odeint, solve_ivp
from Derivative_Equations import *

G = 1.33

# Solving for the derivative at the critical point
# Function to obtain the numerators and denominator at the c.p.
def Crit_derivatives(a, ini_guess):
    def Crit(S, a):
        r, v, rho, Pressure, L = S
        # Br = np.sqrt(4*np.pi*abs(Pressure)/30)
        Br = 0
        Bphi = Br
        eq1, temp = dvdr([r, v, Pressure, rho, L, Br, Bphi, a])
        eq2, den = dPdr([r, v, Pressure, rho, L, Br, Bphi, a])
        eq3, temp = drhodr([r, v, Pressure, rho, L, Br, Bphi, a])
        eq4, temp = dLdr([r, v, Pressure, rho, L, Br, Bphi, a])
        # Use this when using fsolve
        # return [eq1, eq2, eq3, den]
        
        # Use this when using leastsq
        return [eq1, eq2, eq3, eq4, den]

    r_crit, v_crit, rho_crit, P_crit, L_crit = fsolve(Crit, ini_guess, args=(0), maxfev = 20000)
    print(f"r_crit = {r_crit}")
    print(f"v_crit = {v_crit}")
    print(f"P_crit = {P_crit}")
    print(f"rho_crit = {rho_crit}")
    print(f"L_crit = {L_crit}")
    
    Crits = np.array([r_crit, v_crit, P_crit, rho_crit, L_crit])

    Br_crit = np.sqrt(4*np.pi*abs(P_crit)/300)
    #Br_crit = 0
    Bphi_crit = Br_crit

    S = [r_crit, v_crit, P_crit, rho_crit, L_crit, Br_crit, Bphi_crit, 0]
    num1, den1 = dvdr(S)
    num2, den2 = dPdr(S)
    num3, den3 = drhodr(S)
    num4, den4 = dLdr(S)

    print(f"num1 = {num1}, den1 = {den1}")
    print(f"num2 = {num2}, den2 = {den2}")
    print(f"num3 = {num3}, den3 = {den3}")
    print(f"num4 = {num4}, den4 = {den4}")

    # Defining step size
    h = 0.005

    # Calculating derivative at point slightly above critical radius
    def crit_der_up(S, a, Crits):
        dvdr_crit, dPdr_crit, drhodr_crit, dLdr_crit = S
        r_crit, v_crit, P_crit, rho_crit, L_crit = Crits
        r = r_crit + h
        v = v_crit + h*dvdr_crit
        P = P_crit + h*dPdr_crit
        rho = rho_crit + h*drhodr_crit
        L_calc = L_crit + h*dLdr_crit 
        Br = np.sqrt(4*np.pi*abs(P)/300)
        #Br = 0
        Bphi = Br
        eq1 = dvdr([r, v, P, rho, L_calc, Br, Bphi, a])[0]/dvdr([r, v, P, rho, L_calc, Br, Bphi, a])[1] - dvdr_crit
        eq2 = dPdr([r, v, P, rho, L_calc, Br, Bphi, a])[0]/dPdr([r, v, P, rho, L_calc, Br, Bphi, a])[1] - dPdr_crit
        eq3 = drhodr([r, v, P, rho, L_calc, Br, Bphi, a])[0]/drhodr([r, v, P, rho, L_calc, Br, Bphi, a])[1] - drhodr_crit
        eq4 = dLdr([r, v, P, rho, L_calc, Br, Bphi, a])[0]/dLdr([r, v, P, rho, L_calc, Br, Bphi, a])[1] - dLdr_crit
        return [eq1, eq2, eq3, eq4]
    
    # Calculating derivative at point slightly below critical radius
    def crit_der_down(S, a, Crits):
        dvdr_crit, dPdr_crit, drhodr_crit, dLdr_crit = S
        r_crit, v_crit, P_crit, rho_crit, L_crit = Crits
        r = r_crit - h
        v = v_crit - h*dvdr_crit
        Pressure = P_crit - h*dPdr_crit
        rho = rho_crit - h*drhodr_crit
        L_calc = L_crit + h*dLdr_crit     
        Br = np.sqrt(4*np.pi*abs(Pressure)/300)
        #Br = 0
        Bphi = Br
        eq1 = dvdr([r, v, Pressure, rho, L_calc, Br, Bphi, a])[0]/dvdr([r, v, Pressure, rho, L_calc, Br, Bphi, a])[1] - dvdr_crit
        eq2 = dPdr([r, v, Pressure, rho, L_calc, Br, Bphi, a])[0]/dPdr([r, v, Pressure, rho, L_calc, Br, Bphi, a])[1] - dPdr_crit
        eq3 = drhodr([r, v, Pressure, rho, L_calc, Br, Bphi, a])[0]/drhodr([r, v, Pressure, rho, L_calc, Br, Bphi, a])[1] - drhodr_crit
        eq4 = dLdr([r, v, Pressure, rho, L_calc, Br, Bphi, a])[0]/dLdr([r, v, Pressure, rho, L_calc, Br, Bphi, a])[1] - dLdr_crit
        return [eq1, eq2, eq3, eq4]

    dvdr_ini_guess = -0.4
    dPdr_ini_guess = -10
    drhodr_ini_guess = -25
    dLdr_ini_guess = 0.01

    crit_der_ini_guess = [dvdr_ini_guess, dPdr_ini_guess, drhodr_ini_guess, dLdr_ini_guess]

    #Solving for the derivatives above rc
    v_crit_up, P_crit_up, rho_crit_up, L_crit_up = fsolve(crit_der_up, crit_der_ini_guess, args = (a, Crits), maxfev = 2000)
    print(f"Plugging in the derivatives above rc: {crit_der_up([v_crit_up, P_crit_up, rho_crit_up, L_crit_up], a, Crits)}")
    
    #Solving for the derivatives below rc
    v_crit_down, P_crit_down, rho_crit_down, L_crit_down = fsolve(crit_der_down, crit_der_ini_guess, args = (a, Crits), maxfev = 2000)
    print(f"Plugging in the derivatives below rc: {crit_der_down([v_crit_up, P_crit_up, rho_crit_up, L_crit_up], a, Crits)}")

    return [r_crit, v_crit, P_crit, rho_crit, v_crit_up, P_crit_up, rho_crit_up, L_crit_up], [r_crit, v_crit, P_crit, rho_crit, v_crit_down, P_crit_down, rho_crit_down, L_crit_down]

# Function to calculate the derivatives at ANY point
def derivatives(r, S, Lc, a, Crit_vals, B_max):
    r_crit, v_crit, dvdr_critical, dPdr_critical, drhodr_critical, dLdr_critical = Crit_vals
    v, p, rho, l = S
    Br = np.sqrt(4*np.pi*abs(p)/300) * r_crit/r
    Bphi = (np.sqrt(4*np.pi*abs(p)/300)*v_crit - np.sqrt(4*np.pi*abs(p)/30)*Lc/r_crit + Br*l/r) / v
    #Br = B_max * r_crit/r
    #Bphi = (B_max*v_crit - B_max*Lc/r_crit + Br*l/r) / v
    Vars = [r, v, p, rho, l, Br, Bphi, a]
    if (r == r_crit):
        v_grad = dvdr_critical
        P_grad = dPdr_critical
        rho_grad = drhodr_critical
        L_grad = dLdr_critical
    else:
        v_grad = dvdr(Vars)[0] / dvdr(Vars)[1] 
        P_grad = dPdr(Vars)[0] / dPdr(Vars)[1]
        rho_grad = drhodr(Vars)[0] / drhodr(Vars)[1]
        L_grad = dLdr(Vars)[0] / dLdr(Vars)[1]
    return np.array([v_grad, P_grad, rho_grad, L_grad])
    
# Solving the system of ODEs
def Solve_ODE(Lc, a, sgn):
    S_up, S_down = Crit_derivatives(Lc, a, sgn)
    r_crit, v_crit, P_crit, rho_crit, dvdr_critical_up, dPdr_critical_up, drhodr_critical_up, dLdr_critical_up = S_up
    Crits_up = np.array([r_crit, v_crit, dvdr_critical_up, dPdr_critical_up, drhodr_critical_up, dLdr_critical_up])
    
    r_crit, v_crit, P_crit, rho_crit, dvdr_critical_down, dPdr_critical_down, drhodr_critical_down, dLdr_critical_down = S_down
    Crits_down = np.array([r_crit, v_crit, dvdr_critical_down, dPdr_critical_down, drhodr_critical_down, dLdr_critical_down])
    
    Crits = (Crits_down - Crits_up) / 2
    
    #Solving the ODE above rc
    r_solve_above = np.arange(r_crit, 15, 0.01)
    solution_above = solve_ivp(derivatives, t_span = np.array([min(r_solve_above), max(r_solve_above)]), y0 = np.array([v_crit, P_crit, rho_crit, Lc]), method = 'DOP853' , t_eval = r_solve_above, args = (Lc, a, Crits))

    #Solving the ODE Below rc
    r_solve_below = np.arange(2.05, r_crit, 0.001)
    solution_below = solve_ivp(derivatives, t_span = np.array([min(r_solve_below), max(r_solve_below)]), y0 = np.array([v_crit, P_crit, rho_crit, Lc]), method = 'DOP853' , t_eval = r_solve_below, args = (Lc, a, Crits))

    #Reversing the np arrays for the solution below rc
    #r_solve_below = np.flipud(r_solve_below)
    #solution_below = np.flipud(solution_below)

    #Combining the arrays for the solution
    r_solve = np.concatenate((r_solve_below, r_solve_above))
    solution = np.concatenate((solution_below, solution_above))

    #Determining v_r and c_s
    v_r = np.array([])
    P_sol = np.array([])
    rho = np.array([])
    L_sol = np.array([])
    for i in range(len(r_solve)):
        v_r = np.concatenate((v_r, np.array([solution[i][0]])))
        P_sol = np.concatenate((P_sol, np.array([solution[i][1]])))
        rho = np.concatenate((rho, np.array([solution[i][2]])))
        L_sol = np.concatenate((L_sol, np.array([solution[i][3]])))
        
    #For Accretion
    if dvdr_critical_down < 0:
        str = "Accretion"
        x = v_r
        
    else:
        str = "Wind"
        x = np.sqrt(G * P_sol / rho)
    
    c_s = np.sqrt(G * P_sol / rho)
    
    
    #Plotting the final solution
    plt.plot(r_solve, c_s, label = "Speed of sound")
    plt.plot(r_solve, v_r, label = "Speed of gas")
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

    #Plotting the Pressure
    plt.plot(r_solve, P_sol)
    plt.grid()
    plt.xlabel("The distance from center of black hole (in r0 units)")
    plt.ylabel("Pressure of the gas")
    plt.title("Plot of Pressure as a function of r in case of " + str)
    plt.show()

    #Plotting the density
    plt.plot(r_solve, rho)
    plt.grid()
    plt.xlabel("The distance from center of black hole (in r0 units)")
    plt.ylabel("Density of the gas")
    plt.title("Plot of Density as a function of r in case of " + str)
    plt.show()

    #Plotting the angular momentum
    plt.plot(r_solve, L_sol)
    plt.grid()
    plt.xlabel("The distance from center of black hole (in r0 units)")
    plt.ylabel("Angular momentum of the gas")
    plt.title("Plot of Angular momentum as a function of r in case of " + str)
    plt.show()

    final_solution = [v_r, P_sol, rho, L_sol]
    print(f"The critical radius is: {r_crit}")   
    print(f"The critical velocity derivative above is: {dvdr_critical_up}") 
    print(f"The critical velocity derivative below is: {dvdr_critical_down}") 
    print(solution_above[5][0])
    
    return final_solution, r_solve