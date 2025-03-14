from electrokitty import ElectroKitty
from electrokitty.ElectroKitty_parameter_distributions import gaussian_distribution
from electrokitty.ElectroKitty_simulator import python_electrokitty_simulator as ps
import numpy as np
import matplotlib.pyplot as plt

mechanism = "E(1): a*=c*"

mu = 0.15
sigma = 0.1

kins = [[0.5, 10, 0]]

D = []

ic = [
      [0,10**-5],
      []
      ]

iso = [0,0]

cc = [293, 0, 0, 10**-4]

si = [0.0001, 20, 10**-5, 0]

problem1 = ElectroKitty(mechanism)

#problem1.set_data(E, np.zeros(len(E)), t)
problem1.create_simulation(kins, cc, D, iso, si, ic)
problem1.V_potencial(-0.5, 0.5, 0.1, 0., 0, 200)

sim = ps()
sim.give_sim_program(problem1.E_generated, problem1.t)
sim.update_parameters(problem1.mechanism_list, kins, cc, D, iso, si, ic)
i = sim.simulate()
#print(i)
plt.plot(i[1])
plt.show()

#ed, id, td = problem1.simulate()
#problem1.Plot_simulation()

#problem1.Plot_concentration_profile()
#print(problem1.concentration_profile.size)

""" kins = [[0.5, 10, [gaussian_distribution(0, 0.1), 15, -0.35, 0.35, "lin"]]]

D = []

ic = [
      [0,10**-5],
      []
      ]

iso = [0,0]

cc = [293, 10, 6, 10**-4]

si = [0.0001, 20, 10**-5, 0]

problem = ElectroKitty(mechanism)
#problem.set_data(ed, id, td)
#problem.create_simulation(kins, cc, D, iso, si, ic)
#problem.simulate()
#problem.Plot_data()
#plt.figure()
#plt.plot(ed, id)
#plt.plot(ed, problem.current)
#plt.show()
#problem.fit_to_data(N_disp=15, algorithm="CMA-ES", fit_Ru=True, fit_Cdl=True, fit_iso=True)
#problem.sample_parameter_distribution(n_samples=5, N_disp=5)
#problem.print_fitting_parameters()
#print(problem.safety_tuple)
#print(problem.mean_chain) """