from publish import *

# # Computational Mathematics
# ## Problem sheet 3
# ## 2023-09-30
# ## Patrick Farrell (<patrick.farrell@maths.ox.ac.uk>)
# ## Oriel College

import matplotlib.pyplot as plt
import sympy as sp


#
# ***
# ### Question 1 (helping Lagrange).
#


# Generalised coordinates
t = sp.Symbol("t", real=True, nonnegative=True)
theta = sp.Function(r"\theta", real=True)(t)
r = sp.Function("r", real=True, nonnegative=True)(t)

# Physical parameters
m = sp.Symbol("m", real=True, positive=True)
M = sp.Symbol("M", real=True, positive=True)
G = sp.Symbol("G", real=True, positive=True)

# Express Cartesian
x = r * sp.cos(theta)
y = r * sp.sin(theta)

# Kinetic energy
T = sp.simplify(m/2 * (sp.diff(x, t)**2 + sp.diff(y, t)**2))

# Potential energy
V = -G*M*m/r

# Define Lagrangian
L = T - V
render(L, name=r"L(r, \theta)")

# Notice that Lagrangian does not depend on \theta
assert sp.diff(L, theta) == 0
const_from_ELTheta = sp.diff(L, sp.diff(theta, t))
l = sp.Symbol(r"\ell", real=True, positive=True)
fact = sp.solve([const_from_ELTheta - l], [sp.diff(theta, t)])
render(fact)

# Now calculate Euler-Lagrange equation for r
# and substitute in known information about theta
EL = sp.diff(sp.diff(L, sp.diff(r, t)), t) - sp.diff(L, r)
EL = EL.subs(fact)
render(EL, name="0")

#
# ***
# ### Question 2 (helping Pauli).
#

n = sp.Symbol("n", integer=True, positive=True)
l = sp.Symbol(r"\ell", integer=True, nonnegative=True)
a = sp.Symbol("a", real=True, positive=True)
r = sp.Symbol("r", real=True, nonnegative=True)

R = (
      sp.sqrt( (2/(n*a))**3 * sp.factorial(n - l - 1) / (2*n * sp.factorial(n + l)) )
    * sp.exp(-r / (n*a))
    * ((2*r)/(n*a))**l
    * sp.assoc_laguerre(n - l - 1, 2*l + 1, (2*r)/(n*a))
    )

render(R, name="R_{n, \ell}")


def f(n_val, l_val, k):
    """
    Evaluate the integrand for specific values of n, l, k.
    """
    R_nl = R.subs({n: n_val, l: l_val})
    return sp.integrate(R_nl**2 * r**k, (r, 0, sp.oo))


mu_10 = f(1, 0, 3)
sigma_10 = sp.sqrt(f(1, 0, 4) - mu_10**2)
render(mu_10, name=r"\mu_{10}")
render(sigma_10, name=r"\sigma_{10}")


# Main calculation for l = 0
ns = list(range(1, 9))
mus_0 = []
sigmas_0 = []
for n_val in ns:
    mu = f(n_val, 0, 3)
    sigma = sp.sqrt(f(n_val, 0, 4) - mu**2)

    mus_0.append(float(sp.N(mu.subs({a: 1}))))
    sigmas_0.append(float(sp.N(sigma.subs({a: 1}))))


# Main calculation for l = 1
mus_1 = []
sigmas_1 = []
for n_val in ns[1:]:
    mu = f(n_val, 1, 3)
    sigma = sp.sqrt(f(n_val, 1, 4) - mu**2)

    mus_1.append(float(sp.N(mu.subs({a: 1}))))
    sigmas_1.append(float(sp.N(sigma.subs({a: 1}))))


plt.grid()
plt.errorbar(ns, mus_0, fmt='ob', ecolor='k', yerr=sigmas_0, capsize=10, label=r"$\ell = 0$")
plt.errorbar(ns[1:], mus_1, fmt='or', ecolor='k', yerr=sigmas_1, capsize=10, label=r"$\ell = 1$")
plt.xlabel(r"$n$")
plt.ylabel(r"$\mu_{n \ell} \pm \sigma_{n \ell}$ (a)")
plt.title("Distance of the electron to the nucleus of the hydrogen atom")
plt.gcf().set_dpi(300)
plt.xticks(ticks=ns, labels=[str(n_val) for n_val in ns])



#
# ***
#

publish()
