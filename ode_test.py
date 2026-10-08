import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sciint
from electrokitty import ElectroKitty

k1 = 0.5
k2 = 0.

c1 = 10
c2 = 10**-4

# def fun(t, y):
#     return [-k1*y[0]*y[1] + k2*y[1]**2, k1*y[0]*y[1] - k2*y[1]**2]

# # def fun(t, y):
# #     return [-k1*y[0] + k2*y[1], k1*y[0] - k2*y[1]]

# # def analytic_sol(t):
# #     y = (c1 + c2)/((c1/c2 - k2/k1)*np.exp(-k1*(c1 + c2)*t)+1+k2/k1)
# #     x = (c1 + c2)*np.exp(-k1*(c1 + c2)*t)*(c1/c2 - k2/k1)/((c1/c2 - k2/k1)*np.exp(-k1*(c1 + c2)*t)+1+k2/k1)
# #     return x, y

# sol = sciint.solve_ivp(fun, [0, 30], [c1, c2])

# plt.plot(sol.t, sol.y[0], marker = "o")
# plt.plot(sol.t, sol.y[1], marker = "o")

# sols = analytic_sol(sol.t)

# plt.plot(sol.t, sols[0])
# plt.plot(sol.t, sols[1])

sim = ElectroKitty("E(1): a = b")

sim.create_simulation([[0.5, 1000, 0]], [293, 0, 0, 10**-4], [10**-9, 10**-9], [0, 0], [0.001, 25, 10**-5, 0], [[10**-5, 0], [1, 0]])
# sim.C_potential(0, max(sol.t), 1000)
sim.V_potential(0.5, -0.5, 0.05, 0, 0, 1000)

sim.simulate()

sim.Plot_simulation()
print("here")

# plt.plot(sim.t, sim.surface_profile[:, 0])
# plt.plot(sim.t, sim.surface_profile[:, 1])

# sim.Plot_Adsorbed_species()