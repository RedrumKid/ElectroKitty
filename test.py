from electrokitty import ElectroKitty
import numpy as np
import matplotlib.pyplot as plt

mechanism = "C: a* = b* \n E(1):b* = c*"

kin = [[1000, 750],
       [0, 4000, 0.8]]

ic = [[0, 0, 1.2*10**-5], []]

D = []

#iso = [-6, -4, -1, -10]

iso = [[[0, 0, 0], [-3, 2, 0]],
       [[0, -2, 0], [0, 0, -4.5]]]

si = [0.001, 10, 10**-5, 0]

cc = [293, 0, 0.4, 0.125*10**-4]

sim = ElectroKitty(mechanism)

sim.V_potential(0.45, 0.875, 0.05, 0, 0, 1000)

sim.create_simulation(kin, cc, D, iso, si, ic)

e, i, t = sim.simulate()
sim.Plot_Adsorbed_species()
plt.plot(e, i)
plt.show()