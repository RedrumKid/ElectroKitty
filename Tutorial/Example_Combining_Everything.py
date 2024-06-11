# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 16:17:42 2024

@author: ozbejv
"""

import numpy as np
import matplotlib.pyplot as plt
from ElectroKitty_main import simulator_Main_loop, V_potencial, uniform

# Constants used in simulation
# All constants are given in base SI units (V, A, m, s, mol)

F=96485 #As/mol
R=8.314 #J/mol/K
T=293 #K

# Here I combine all the previous steps to show how to simulate species
# that can simultaniouslly adsorb and react in the bulk
# I have also simulate repeling forces between the adsorbed species, by setting 
# the isotherm constants as negative, a positive number would mean attraction

# Number of points in the spatial direction
nx=20

# Mechanism string written in full
mechanism="E(2):a=b \n C: b=c \n C: b+*-b* \n E(1):b*=d*"
# Number of dissolved species, must be updated so to not brick the simulation

# Constants given for adsorption and reaction on the surface
# One constant for ireversible, two for reversible
# units are as the apply, some examples are given
kin_const=[
    [0.5,10**2,0], # /, m/s, V
    [10,1], # /s
    [10],
    [0.5,100,-0.15]
    ]

# The diffusion constants for the dissolved species
# This can be any list containing the constants, must be of lenghth n_spec
D=3*[10**-9]

viscosity=10**-5 #m^2/s
rot_freq=0 #Hz

# Constants describing the cell: Temperature, Ru, Cdl, A
cell_c=[T, 0, 1*10**-4, 10**-4] # K, Ohm , F/m^2, m^2

# Constants describing the x direction: dx/xmax, nx
si=[0.01/36, nx, viscosity, rot_freq]


# Initial condition: first list initital surface occupations, second a list of initial condition functions, reccomend use of the uniform function 
# uniform function takes in the value of the initial concentration and give the correct way for the simulator to understand the given list
spec_info=[
    [10**-5,0,0], #[mol/m^2]
    [uniform(1),uniform(0),uniform(0)] #[mol/m^3]
    ]

isotherm=[-0.1,-0.1,-0.1]

# Constants used in describing the potential program
# In this case I use the CV programm, other are to follow, a custom one can be made
# by generating a the time domain and the electrical potential as numpy arrays

Ei=0.3  #Initial potential, also the potential to wich the CV returns [V]
Ef=-0.3 #Final potential, we cycle to this potential [V]
v=0.1   #scan speed [V/s]
freq=9  #In case of ACV this is the frequency of a sine wave [Hz]
amp=0.0  #In case of ACV this is the amplitude of a sine wave [V]
nt=500 #number of time points

# Generate the potential program
E,t=V_potencial(Ei,Ef,v,amp,freq,nt,F/R/T) # V, s

# Run the simulation
E_corr,i,t=simulator_Main_loop(mechanism, [kin_const, cell_c, D, isotherm], si, t, spec_info, E) #V, A, sA


# A simple plot of the calculated current
plt.figure("Example plot", figsize=(9,5))
plt.title("Example CV")
plt.plot(E,i, linewidth=4)
plt.ylabel("i [A]")
plt.xlabel("E [V]")

plt.show()
