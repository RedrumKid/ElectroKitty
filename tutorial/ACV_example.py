# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 15:12:03 2024

@author: ozbejv
"""

from electrokitty import ElectroKitty
import numpy as np
import matplotlib.pyplot as plt

# Simmilairly as before we start by creating a mechanism and passing the data to the class

# Proceed as in the EC example

mechanism = "E(1): a=b \n C: b>c"

problem = ElectroKitty(mechanism)

kinetic_cons = [
    [0.51, 0.1, 0.0],
    [0.1]
    ]

D = 3*[10**-9]

initial_concentrations = [
    [],
    [1,0,0]
]

isotherm_constants = []

solution_viscosity = 10**-5
rotational_frequency = 0

cell_constants = [298, 00, 0*10**-6, 10**-4]

spatial_information = [0.01/36, 20, solution_viscosity, rotational_frequency]

Ei = 0.5
Ef = -0.5
v = 0.05
amp = 0
freq = 9
nt = 500

data = np.loadtxt("test_ACV.txt")
i_data = data[:,1]

problem.set_data(data[:,0], i_data, data[:,2])

problem.create_simulation(kinetic_cons, cell_constants, 
                          D, isotherm_constants, 
                          spatial_information, initial_concentrations)

# problem.simulate()
problem.Plot_data()
# problem.Plot_simulation()    

# Here we define the number of harmonics we want to plot
# and the a list of widths in the FT of the data to filter the individual harmonics

N_harmonics = 5
w = 0.1*freq*np.ones(N_harmonics+1)

# calling the function to analyse data
problem.FFT_analyze_data(freq, N_harmonics, w)
# plotting the data
problem.Harmonic_plots(plot_data = True)

# the same here only that we plot the simulated data
problem.FFT_analyze_sim(freq, N_harmonics, w)
problem.Harmonic_plots(plot_sim = True)