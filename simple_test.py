from electrokitty import ElectroKitty
from electrokitty.ElectroKitty_parameter_distributions import gaussian_distribution
import numpy as np
import matplotlib.pyplot as plt

mechanism = "E(1): a*=c* \n E(1): c*=b*"

mu = 0.15
sigma = 0.1

kin = [[0.5, 10, [gaussian_distribution(mu, sigma), 15, -3.5*sigma+mu, 3.5*sigma+mu, "lin"]], 
        [0.5, 10, -0.15]]

D = []

ic = [
      [0, 0, 10**-5],
      []
      ]

iso = [0,0,0]

cc = [293, 0, 0, 10**-4]

si = [0.0001, 20, 10**-5, 0]

problem = ElectroKitty(mechanism)
problem.V_potencial(-0.5, 0.5, 0.1, 0., 0, 500)
# problem.set_data(E, np.zeros(len(E)), t)
problem.create_simulation(kin, cc, D, iso, si, ic)

ed, id, td = problem.simulate()
plt.figure()
plt.plot(ed, id)
#problem.Plot_simulation() """

""" print(problem.simulator.simulate_with_dispersion)
print(problem.simulator.dispersed_cell_const)
print(problem.simulator.dispersed_isotherm)
print(problem.simulator.dispersed_kin)
print(problem.simulator.dispersed_species_information) """


mechanism = "E(1): a*=c* \n E(1): c*=b*"

mu = 50
sigma = 1

kin = [[0.5, 10, 0.15], 
        [0.5, 10, -0.15]]

D = []

ic = [
      [0, 0, 10**-5],
      []
      ]

iso = [0,0,0]

cc = [293, 0, 1.65, 10**-4]

si = [0.0001, 20, 10**-5, 0]

problem = ElectroKitty(mechanism)
problem.V_potencial(-0.5, 0.5, 0.1, 0., 0, 500)
# problem.set_data(E, np.zeros(len(E)), t)
problem.create_simulation(kin, cc, D, iso, si, ic)

e, i, t = problem.simulate()
plt.plot(e, i)
#problem.Plot_simulation()

cc = [293, 0, 10**0, 2*10**-4]

si = [0.0001, 20, 10**-5, 0]

plt.show()