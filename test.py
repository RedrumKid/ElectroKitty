from electrokitty import ElectroKitty
import numpy as np
import matplotlib.pyplot as plt

mechanism = "E(1): a* = b* \n E(1): b* = c*"

kin = [[0.5, 1000, 0.4],
       [0.5, 1000, -0.4]]

ic = [[10**-5, 0, 0], []]

D = []

iso = [[[-6, -6, 0], [-6, -6, 0]],
       [[0., -1, -1], [0., -1, -1]]]

si = [0.001, 10, 10**-5, 0]

cc = [293, 0, 0, 10**-4]

sim = ElectroKitty(mechanism)

sim.V_potential(0.75, -0.75, 0.1, 0, 0, 1000)

sim.create_simulation(kin, cc, D, iso, si, ic)

e, i, t = sim.simulate()

print(i)

plt.plot(e, i)
plt.show()