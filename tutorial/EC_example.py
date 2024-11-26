# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 15:12:03 2024

@author: ozbejv
"""

from electrokitty import ElectroKitty
import numpy as np
import matplotlib.pyplot as plt

# Simmilairly as before we start by creating a mechanism and passing the data to the class

# Proceed as before

# in the mechanism we separate the steps with a \n, and each step must be declared as either E or C
# The =  and > or < are showing the direction of the reaction 
# Unlike the Surface confined example we leave out the * at the end to declare the species as dissolved
mechanism = "E(1): a=b \n C: b>c"

problem = ElectroKitty(mechanism)

kinetic_cons = [
    [0.51, 0.1, 0.0],
    [0.1]
    ]

D = 3*[10**-9]

# Here we start filling the second list and leave the first one empty 

initial_concentrations = [
    [],
    [1,0,0]
]

isotherm_constants = []

solution_viscosity = 10**-5
rotational_frequency = 0

# Temperature, uncompensated resistance, double-layer capacitance, area
# K, Ohm, F, m^2
cell_constants = [298, 00, 0*10**-6, 10**-4]

# initial fraction, number of points
spatial_information = [0.01/36, 20, solution_viscosity, rotational_frequency]

Ei = 0.5
Ef = -0.5
v = 0.05
amp = 0
freq = 0
nt = 500

data = np.loadtxt("test.txt")
i_data = data[:,1]
noise = np.random.normal(0, 0.01*max(i_data), len(data[:,1]))
i_data += noise
problem.set_data(data[:,0], i_data, data[:,2])

problem.create_simulation(kinetic_cons, cell_constants, 
                          D, isotherm_constants, 
                          spatial_information, initial_concentrations)

problem.simulate()
problem.Plot_data()
problem.Plot_simulation()    
    