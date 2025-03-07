from electrokitty import ElectroKitty
from electrokitty.ElectroKitty_parameter_distributions import gaussian_distribution
import numpy as np
import matplotlib.pyplot as plt

""" mechanism = "E(1): a*=c*"

mu = 0.15
sigma = 0.1

kins = [[0.5, 10, 0]]

D = []

ic = [
      [0, 10**-5],
      []
      ]

iso = [0,0]

cc = [293, 0, 0, 10**-4]

si = [0.0001, 20, 10**-5, 0]

problem1 = ElectroKitty(mechanism)

# problem.set_data(E, np.zeros(len(E)), t)
problem1.create_simulation(kins, cc, D, iso, si, ic)
problem1.V_potencial(-0.5, 0.5, 0.1, 0., 0, 500)

ed, id, td = problem1.simulate()
problem1.save("test")
plt.figure()
plt.plot(ed, id) """

problem = ElectroKitty()
problem.load("5_index_good_fit.ek")
problem.Plot_data()
print(problem.kin)
plt.show()