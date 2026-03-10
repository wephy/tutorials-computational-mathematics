# # Computational Mathematics
# ## Problem sheet 3
# ## 2026-01-26
# ## Alex Hanson (<worc6973@ox.ac.uk>)
# ## Worcester College

#
# ***
# ### Question 1.
#


import sympy as sp

t = sp.symbols('t', real = True)
m, M, G = sp.symbols('m M G', positive = True)

r  = sp.Function('r')(t)
theta = sp.Function('theta')(t)

rdot = sp.diff(r, t)
thetadot = sp.diff(theta, t)

T = sp.Rational(1,2)*m*(rdot**2 + r**2*thetadot**2)
V = -G*M*m/r

L = T - V
sp.simplify(L)
print(L)

dLdthetadot = sp.diff(L,thetadot)
sp.simplify(dLdthetadot)

#Since L does not depend on theta, dL/dtheta = 0

#
angularmomentum = sp.symbols('ang_momentum', positive=True)

theta_dot_expression = sp.solve(sp.Eq(m*r**2*thetadot, angularmomentum), thetadot)[0]

print(theta_dot_expression)

dLdrdot_t = sp.diff(sp.diff(L, rdot), t) - sp.diff(L, r)
#so dLdrdot_t = 0
sp.simplify(dLdrdot_t)

dLdrdot_t_sub = dLdrdot_t.subs(thetadot, theta_dot_expression)
sp.simplify(dLdrdot_t_sub)

rdoubledot = sp.diff(rdot, t)
solution = sp.solve(dLdrdot_t_sub, rdoubledot)[0]
sp.simplify(solution)

print(solution)

#Outputs:
#G*M*m/r(t) + m*(r(t)**2*Derivative(theta(t), t)**2 + Derivative(r(t), t)**2)/2
#ang_momentum/(m*r(t)**2)
#-G*M/r(t)**2 + ang_momentum**2/(m**2*r(t)**3)