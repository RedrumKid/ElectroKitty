import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sciint
from electrokitty import ElectroKitty

k1 = 0.5
k2 = 0.

c1 = 10
c2 = 1

def fun(t, y):
    return [-k1*y[0] + k2*y[1]**2, 2*k1*y[0] - 2*k2*y[1]**2]

# def fun(t, y):
#     return [-k1*y[0] + k2*y[1], k1*y[0] - k2*y[1]]

def analytic_sol(t):
    y = (c1 + c2)/((c1/c2 - k2/k1)*np.exp(-k1*(c1 + c2)*t)+1+k2/k1)
    x = (c1 + c2)*np.exp(-k1*(c1 + c2)*t)*(c1/c2 - k2/k1)/((c1/c2 - k2/k1)*np.exp(-k1*(c1 + c2)*t)+1+k2/k1)
    return x, y

sol = sciint.solve_ivp(fun, [0, 100], [c1, c2])

plt.plot(sol.t, sol.y[0], marker = "o")
plt.plot(sol.t, sol.y[1], marker = "o")

# sols = analytic_sol(sol.t)

# plt.plot(sol.t, sols[0])
# plt.plot(sol.t, sols[1])

sim = ElectroKitty("C: a* = b* + b*")

sim.create_simulation([[k1, k2]], [293, 0, 0, 1], [], [0, 0], [0.1, 20, 10**-5, 0], [[c1, c2], []])
sim.C_potential(0, max(sol.t), 1000)

sim.simulate()

plt.plot(sim.t, sim.surface_profile[:, 0])
plt.plot(sim.t, sim.surface_profile[:, 1])

# sim.Plot_Adsorbed_species()