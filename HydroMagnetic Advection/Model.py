import sympy as smp

#Solving for the equations of motion of the system
r, F, h, dhdr, rho, drhodr, v, dvdr = smp.symbols('r, F, h, dhdr, rho, drhodr, v, dvdr')
P, dPdr, L, dLdr, Br, dBrdr, BPhi, dBPhidr = smp.symbols('P, dPdr, L, dLdr, Br, dBrdr, BPhi, dBPhidr')
G, alpha, fvis, fm, dFdr = smp.symbols('G, alpha, fvis, fm, dFdr')

##Set of equations that are to be solved
eqn1 = Br + r*dBrdr
eqn2 = dvdr*BPhi + v*dBPhidr + Br*L/r**2 - dBrdr*L/r - Br*dLdr/r
eqn3 = rho*h*v + r*drhodr*h*v + r*rho*dhdr*v + r*rho*h*dvdr
eqn4 = v*dvdr + dPdr/rho - L**2/r**3 - (Br*dBrdr - BPhi**2/r)/(4*smp.pi*rho) + F
eqn5 = (2*alpha*(P + rho*v**2) + r*alpha*(dPdr + drhodr*v**2 + 2*v*rho*dvdr)) / rho + (Br*dBrdr + Br*BPhi/r) / (4*smp.pi*rho) - v*dLdr
eqn6 = v*(dPdr - G*P*drhodr/rho)/(G - 1) - alpha*fvis*(P + rho*v**2)*L/r**2 + 3*fm*(Br**2 + BPhi**2) * v/(16*smp.pi*r) 
eqn7 = h*(1/r + (dPdr)/(P) - dFdr/F - drhodr/rho) / 2 - dhdr

##Solving the above equations to obtain the expressions for the derivatives
eqns = [eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7]
sol = smp.linsolve(eqns, [dvdr, dPdr, drhodr, dLdr, dBrdr, dBPhidr, dhdr])

print(sol)
