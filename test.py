import numpy as np
import matplotlib.pyplot as plt
from electrokitty import ElectroKitty

mechanism = "E(2): c = d"

kinetic_constants = [[0.5, 10, 0.]] # A list contaning a list of constants in the order: alpha [/], k0 [m/s], E0 [V]

Diffusion_constants=[10**-9, 10**-9, 10**-9, 10**-9] # The diffusion constats the first is for a, the second for b [m^2/s]

initial_condition = [[], [1, 0]] # The first list for adsorbed species in mol/m^2, the second for dissolved [mol/m^3] 

spatial_information = [0.001, 20, 10**-5, 0] # constant for grid generation and calculating solution velocities for the RDE electrode
# the fraction of the initial step, the number of points for the simulation, solution viscosity [m^2/s] and rotational frequency [Hz]

cell_constants = [293, 0, 0, 0.283*10**-4] # the cell constants: temperature [K], uncompensated resistance [ohm], double-layer capacitance [F/m^2], electrode area[m^2]

iso = [] # we leave it alone

sim = ElectroKitty(mechanism)

sim.create_simulation(kinetic_constants, cell_constants, Diffusion_constants, iso, spatial_information, initial_condition)

sim.V_potential(0.5, -0.5, 0.05, 0, 0, 1000)

e, i, t = sim.simulate()

plt.plot(e, i)