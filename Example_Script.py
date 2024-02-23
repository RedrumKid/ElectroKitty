# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 16:17:42 2024

@author: ozbejv
"""

import numpy as np
import matplotlib.pyplot as plt
from ElectroKitty_main import simulator_Main_loop, V_potencial, uniform, FFT_analysis, Harmonic_plots

# Constants used in simulation
# All constants are given in base SI units (V, A, m, s, mol)

F=96485 #As/mol
R=8.314 #J/mol/K
T=293 #K

# Number of points in the spatial direction
nx=20

# Mechanism string written in full
mechanism="E:a=b "
# Number of dissolved species, must be updated so to not brick the simulation
num_bulk_spec=1

# Constants given for adsorption and reaction on the surface
# One constant for ireversible, two for reversible
Ads_c=[
        # [10]
        ]
# Constants for the bulk reaction
B_c=[
        # [10,1]
      ]
# Electrochemical constants: alpha, k0, E0'
# Electrochemical constant always take 3, regerdles of reversibility
ec_c=[
      [0.5,10**2,0], # none, m/s, V 
      # [0.5,100,-0.15]
      ]

# The diffusion constants for the dissolved species
# This can be any list containing the constants, must be of lenghth n_spec
D=10**-9*np.ones(num_bulk_spec)

# Constants describing the cell: Temperature, Ru, Cdl, A
cell_c=[T, 30, 1*10**-4, 10**-4] # K, Ohm , F/m^2, m^2

# Constants describing the x direction: dx/xmax, nx
si=[0.1/36, nx]


# Initial condition: first list initital surface occupations, second a list of initial condition functions, reccomend use of the uniform function 
# uniform function takes in the value of the initial concentration and give the correct way for the simulator to understand the given list
spec_info=[
    [], #[mol/m^2]
    [uniform(1)] #[mol/m^3]
    ]

# Constants used in describing the potential program
# In this case I use the CV programm, other are to follow, a custom one can be made
# by generating a the time domain and the electrical potential as numpy arrays

Ei=0.3  #Initial potential, also the potential to wich the CV returns [V]
Ef=-0.3 #Final potential, we cycle to this potential [V]
v=0.1   #scan speed [V/s]
freq=9  #In case of ACV this is the frequency of a sine wave [Hz]
amp=0.  #In case of ACV this is the amplitude of a sine wave [V]
nt=1000 #number of time points

# Generate the potential program
E,t=V_potencial(Ei,Ef,v,amp,freq,nt,F/R/T) # V, s

# Run the simulation
E,i,t=simulator_Main_loop(mechanism, [Ads_c, B_c, ec_c, cell_c, D], si, t, spec_info, E) #V, A, s

# A simple plot of the calculated current
plt.figure("Example plot", figsize=(9,5))
plt.title("Example CV")
plt.plot(E,i)
plt.ylabel("i [A]")
plt.xlabel("E [V]")
plt.show()
