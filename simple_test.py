from electrokitty import ElectroKitty
from electrokitty.ElectroKitty_parameter_distributions import gaussian_distribution
import numpy as np
import matplotlib.pyplot as plt

mechanism = "E(1): a*=c*"

mu = 0.15
sigma = 0.1

kins = [[0.5, 10, [gaussian_distribution(0, 0.1), 15, -0.35, 0.3, "lin"]]]

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
plt.figure()
plt.plot(ed, id)

kins = [[0.5, 10, [gaussian_distribution(0, 0.1), 15, -0.35, 0.3, "lin"]]]

D = []

ic = [
      [0, 10**-5],
      []
      ]

iso = [0,0]

cc = [293, 0, 0, 10**-4]

si = [0.0001, 20, 10**-5, 0]

problem = ElectroKitty(mechanism)
#problem.V_potencial(-0.5, 0.5, 0.1, 0., 0, 500)
# problem.set_data(E, np.zeros(len(E)), t)
problem.set_data(ed, id, td)
problem.create_simulation(kins, cc, D, iso, si, ic)


#e, i, t = problem.simulate()
#problem.fit_to_data()
problem.sample_parameter_distribution(n_samples=3, fit_gamamax=True)
#plt.plot(problem.E_generated, problem.current)
#problem.Plot_simulation()
problem.Plot_MCMC_mean_chain()

plt.show()