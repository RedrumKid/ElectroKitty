from electrokitty import ElectroKitty
import numpy as np
import matplotlib.pyplot as plt

mechanism = "E(1): c*=b* \n E(1): b*=a*"

kin = [[0.5, 1, -0.05], [0.5, 10, 0.0]]

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

problem.simulate()
problem.Plot_simulation()