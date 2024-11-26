# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 14:54:46 2024

@author: ozbejv
"""

import numpy as np
import matplotlib.pyplot


# we import the base class as that is all that we will need in our tutorial
from electrokitty import ElectroKitty

# we define our mechanism in this case it is a simple 1e surface confined reaction
mechanism = "E(1): a*=b*"

#calling the class and passing it our mechanism
problem = ElectroKitty(mechanism)

# define our kinetic constants, noticew that it is a list of lists, here it is only 1 
# of lenghth 3. We would need either 1 or 2 in case of a chemical reaction
kinetic_cons = [
    [0.51, 0.1, 0.01],
    ]

# aDefining the diffusion constants of our species, as our species do not diffuse we leave the list empty
D=[]

# Here we define out initial conditions. First is the surfacre bound species, 
# second list contains the initial concentrations

# mol/m2
# mol/m3 
# both can be in whatever unit, as long as those units are consistent throughut the model setup
initial_concentrations = [
    [10**-5,0],
    []
]

# Defining the isotherm constants for each species, dimensionless
# must be defined, 0 is the Langmuir case, any other assumes Frumkin
isotherm_constants = [0,0]

# Defining some cell parameters, the simulator needs both
# Viscosity must not be set ot 
solution_viscosity = 10**-5 # m2/s
rotational_frequency = 0 # Hz

# Temperature, uncompensated resistance, double-layer capacitance, area
# K, Ohm, F, m^2
cell_constants = [298, 00, 0*10**-6, 10**-4]

# initial fraction, number of points
spatial_information = [0.01/36, 20, solution_viscosity, rotational_frequency]

# parameters for creating a CV or ACV potential signal
Ei = 0.5 # initial potential 
Ef = -0.5   # vertex potential
v = 0.05    # scan rate
amp = 0 # sine amplitude
freq = 0    # sine frequency
nt = 500    # number of time points

# importing data from a .txt file
data = np.loadtxt("test_SC.txt")
i_data = data[:,1]
# adding some noise
noise = np.random.normal(0, 0.01*max(i_data), len(data[:,1]))
i_data += noise
# Giving electrokitty the data to use
# The potential for simulation is the one from data
problem.set_data(data[:,0], i_data, data[:,2])
# can use this function to use for a simple simulation

# problem.V_potencial(Ei, Ef, v, amp, freq, nt)

# This function is used for CA measuerments
# Hold potential, maximum time
# problem.C_potential(0, 200, nt)

# setting up the simulation, simply passing the lists to the class
problem.create_simulation(kinetic_cons, cell_constants, 
                          D, isotherm_constants, 
                          spatial_information, initial_concentrations)
# invoke the simulator
problem.simulate()

# plotting
problem.Plot_data()
problem.Plot_simulation()

# uncomment to fit to the imported data, will not work if data not available
# problem.fit_to_data(fit_Ru=True, fit_Cdl=True, fit_iso=True,algorithm="CMA-ES", tolf=10**-11, tolx=10**-11)
# problem.sample_parameter_distribution(n_samples=1000,fit_gamamax=False, fit_Cdl=False, fit_Ru=False)

# plotting the optimised soultion
problem.Plot_simulation()

# plots for the MCMC chain

# problem.Plot_MCMC_parameter_dist()
# problem.Plot_MCMC_histogram()
# problem.Plot_MCMC_mean_chain()

# saving progress as a .ek file, that can be imported via the load function
problem.save("tutorial_SC")