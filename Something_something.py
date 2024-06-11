# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 16:17:42 2024

@author: ozbejv
"""

import numpy as np
import matplotlib.pyplot as plt
from ElectroKitty_main import simulator_Main_loop, V_potencial, uniform, FFT_analysis, Harmonic_plots, Parse_mechanism
import time

# Constants used in simulation
# All constants are given in base SI units (V, A, m, s, mol)

F=96485 #As/mol
R=8.314 #J/mol/K
T=298 #K

# Number of points in the spatial direction
nx=30

# Mechanism string written in full
mechanism="E(1):H_aq+*=H* \n C:H*+H*-*+*"
# Number of dissolved species, must be updated so to not brick the simulation
num_bulk_spec=1

# Constants given for adsorption and reaction on the surface
# One constant for ireversible, two for reversible
# units are as the apply, some examples are given
kin_const=[
    [0.5,10**3,-0.2], # /, m/s, V
    [10**12], # /s
    # [10],
    # [0.5,100,-0.15]
    ]

# The diffusion constants for the dissolved species
# This can be any list containing the constants, must be of lenghth n_spec
D=10**-8*np.ones(num_bulk_spec)

viscosity=10**-5 #m^2/s
rot_freq=60 #Hz

# Constants describing the cell: Temperature, Ru, Cdl, A
cell_c=[T, 0, 0*1*10**-4, 10**-4] # K, Ohm , F/m^2, m^2

# Constants describing the x direction: dx/xmax, nx
si=[0.1/36, nx, viscosity, rot_freq]


# Initial condition: first list initital surface occupations, second a list of initial condition functions, reccomend use of the uniform function 
# uniform function takes in the value of the initial concentration and give the correct way for the simulator to understand the given list
spec_info=[
    [0.75*10**-7, 0], #[mol/m^2]
    [uniform(100)] #[mol/m^3]
    ]

iso=[0,0,0]

# Constants used in describing the potential program
# In this case I use the CV programm, other are to follow, a custom one can be made
# by generating a the time domain and the electrical potential as numpy arrays

Ei=0.1  #Initial potential, also the potential to wich the CV returns [V]
Ef=-0.1 #Final potential, we cycle to this potential [V]
v=0.04   #scan speed [V/s]
freq=1  #In case of ACV this is the frequency of a sine wave [Hz]
amp=0.  #In case of ACV this is the amplitude of a sine wave [V]
nt=200 #number of time points

# Generate the potential program
E,t=V_potencial(Ei,Ef,v,amp,freq,nt,F/R/T) # V, s
# Edc,tdc=V_potencial(Ei,Ef,v,0,freq,nt,F/R/T)
# time1=time.time()
# Run the simulation
E1,i,t=simulator_Main_loop(mechanism, [kin_const, cell_c, D, iso], si, t, spec_info, E) #V, A, s
# # spec, index, types=Parse_mechanism(mechanism)
# # print(time.time()-time1)
# # A simple plot of the calculated current
# # plt.figure("Example plot", figsize=(9,5))

# plt.figure(1)
# plt.plot(E,i)
# plt.plot(E1,i)
# plt.figure(2)
plt.title("Tafel")
plt.plot(np.abs(i[1:]),E1[1:])

