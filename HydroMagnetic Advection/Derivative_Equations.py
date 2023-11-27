import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import sympy as smp

#Defining the variables
r_sym, a_sym = smp.symbols('r, a')
num = (r_sym**2 - 2*a_sym*smp.sqrt(r_sym) + a_sym**2)**2
den = r_sym**3 * (smp.sqrt(r_sym)*(r_sym - 2) + a_sym)**2
force_sym = num/den

#Creating the force function
Force = smp.lambdify((r_sym, a_sym), force_sym)
dFdr = smp.lambdify((r_sym, a_sym), smp.diff(force_sym, r_sym))

#Defining the parameters
G = 1.4
fm = 0.5
alpha = 0
fvis = 0.5
# Lc = 2.5

# The velocity gradient equation
# Function returns numerator and denominator as an array
def dvdr(S):
   r, v, P, rho, L, Br, Bphi, a = S
   num1 = dFdr(r, a)/(2*Force(r, a)) - 3/(2*r) - rho * (1 + 1/G) * (L**2/r**3 - Force(r, a) - (Br**2 + Bphi**2)/(4*np.pi*rho*r)) / (2*P)
   num2 = (G - 1) * (alpha*fvis*(P + rho*v**2)*L/r**2 + 3*fm*(Br**2 + Bphi**2)*v/(16*np.pi*r))
   den = 1/v - v*rho*(1 + 1/G) / (2*P)
   return [- (num1 + num2) * 2*P*G / (rho*v) , - den * 2*P*G / (rho*v)]

# The Pressure gradient equation
# Function returns numerator and denominator as an array
def dPdr(S):
   r, v, P, rho, L, Br, Bphi, a = S
   num1 = G*P*(Br**2 + Bphi**2) / (2*np.pi*rho*v**2*r) + (G - 1) * (L*alpha*fvis*(P + rho*v**2) / (r**2 * v) - 3*fm*(Br**2 + Bphi**2) / (16*np.pi*r))
   num2 = 2*G*P*(Force(r, a) - L**2/r**3) / v**2 - 3*G*P/r + G*P*dFdr(r, a) / Force(r, a)
   den = (G + 1) - 2*G*P/(rho * v**2)
   return [(num1 + num2) , den]

# The Density gradient equation
# Function returns numerator and denominator as an array
def drhodr(S):
   r, v, P, rho, L, Br, Bphi, a = S
   num1 = (G - 1) * (3*fm*(Br**2 + Bphi**2)*rho*(rho * v**2/P - 2)/(16*np.pi*r) + fvis*L*alpha* rho**2 *v*(1 - rho* v**2 / P) / r**2 + 2*L*alpha*fvis*P*rho / (G* r**2))
   num2 = 2* rho**2 * (Force(r, a) - L**2/r**3) - 3* rho**2 * v**2 / r + dFdr(r, a)* rho**2 * v**2 / Force(r, a)
   den = (G + 1) * rho * v**2 - 2*G*P
   return [(num1 + num2)/(rho* v**2) , den/(rho* v**2)]

# The Angular Momentum gradient equation
# Function returns numerator and denominator as an array
def dLdr(S):
   r, v, P, rho, L, Br, Bphi, a = S
   num1 = G*(Bphi**2 + Br**2)*alpha*(P - rho* v**2) / (2*np.pi*r*rho*v) - (G - 1)*3*fm*(Br**2 + Bphi**2)*alpha*v*(1 - rho* v**2 / P) / (16*np.pi*r)
   num2 = G*(Br**2 - Bphi*Br)*P/(2*np.pi* r**2 * rho*v) + (G + 1) * (Bphi*Br + Br**2)*v / (4*np.pi* r**2)
   num3 = 2*G*P*alpha * (Force(r, a) + L**2 * (rho * v**2 / P - 1) / r**3 + v**2 * (rho* v**2 / P + 1/G) / r)
   num4 = (G - 1) * L*alpha*fvis* (P**2 - rho**2 * v**4) / (P* r**2) - 4*G*P*alpha*v/r
   num5 = alpha*v * (rho* v**2 / r + G*P*dFdr(r, a) / Force(r, a) - dFdr(r, a) * rho * v**2 / Force(r, a))
   den = (G + 1) * rho * v**2 / r - 2*G*P/r
   return [(num1 + num2 + num3 + num4 + num5) * r / (rho* v**2) , den * r / (rho* v**2)]