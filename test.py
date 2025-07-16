from electrokitty import ElectroKitty
import matplotlib.pyplot as plt

mechanism = "E(1): a*=c* + c* +c* "

sim = ElectroKitty(mechanism)
sim.create_simulation([[0.5, 100, 0]], [293, 0, 0, 10**-4], [], [0, 0], [0.1, 10, 10, 0], [[1, 0], []])
sim.V_potential(0.5, -0.5, 0.05, 0, 0, 1000)

sim.simulate()

plt.plot(sim.E_generated, sim.surface_profile)
plt.show()