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
T=298 #K

# Number of points in the spatial direction
nx=20

# Mechanism string written in full
# Each step is declared either with E or C, wheather the step is an electrochemical or chemical
# If declared as E, in parenthesis must be an integer, to declare the number of electrons that are used in the reaction 
# eg. (1) for 1 electron
# This is followed by :, and then the simbols for the species. If a species is declared with a * at the end it is adsorbed on the surface
# To denote the arrows use = or -, for a reversible or irreversible reaction. 
mechanism="E(1):a=b"

# Constants given for adsorption and reaction on the surface
# One constant for ireversible, two for reversible
# units are as the apply, some examples are given
alpha=0.5 # /
k0=10**2 #m/s
E0=0 # V

kin_const=[
    [alpha,k0,E0],
    ]

# The diffusion constants for the dissolved species
# This can be any list containing the constants, must have a lenghth equal to the number of dissolved species

D=[10**-9, 10**-9] #m^2/s

# To simulate rotation declare the soulution viscosity and rotation rate
# The rotation rate must be given is Hz instead of the usual RPM
viscosity=10**-5 #m^2/s
rot_freq=0 #Hz

# Constants describing the cell: Temperature, Ru, Cdl, A
Ru=10 #Ohm
Cdl=10**-4 #F
A=10**-4 #m^2
cell_c=[T, Ru, Cdl, A]

# Constants describing the x direction: dx/xmax, nx, solution viscosity and rotation rate of RDE
si=[0.1/36, nx, viscosity, rot_freq]

# Initial condition: first list initital surface concentrations, second a list of initial condition functions for dissolved species,
# reccomend use of the uniform function 
# uniform function takes in the value of the initial concentration and give the correct way for the simulator to understand the given list
c0=1

spec_info=[
    [], #[mol/m^2]
    [uniform(c0),uniform(0)] #[mol/m^3]
    ]

#This is a list for declaring isotherm constants. They will not be used here
isotherm=[] # /

# Constants used in describing the potential program
# In this case I use the CV programm, others are to follow, a custom one can be made
# by generating a the time domain and the electrical potential as numpy arrays

Ei=0.3  #Initial potential, also the potential to wich the CV returns [V]
Ef=-0.3 #Final potential, we cycle to this potential [V]
v=0.1   #scan speed [V/s]
freq=9  #In case of ACV this is the frequency of a sine wave [Hz]
amp=0.08  #In case of ACV this is the amplitude of a sine wave [V]
nt=2**13 #number of time points

# Generate the potential program
# Takes in The initial, final potential, scan rate,
# the amplitude and frequency of the superimposed sin wave, if doing ACV
# The fianl two arguments are the number of points and the electrochemical constant f [1/V]
E,t=V_potencial(Ei,Ef,v,amp,freq,nt,F/R/T) # V, s
# Run the simulation
E_corr,i,t=simulator_Main_loop(mechanism, [kin_const, cell_c, D, isotherm], si, t, spec_info, E) #V, A, s

# A simple plot of the calculated current
plt.figure("Example plot", figsize=(9,5))
plt.rc("font", size=14)
plt.title("Example CV")
plt.plot(E,i, linewidth=4)
plt.ylabel("i [A]")
plt.xlabel("E [V]")

plt.show()


# Here I have used inbuilt functions to separate the harmonics and plot them
# The FFT analysis function requiers the array of the full experiment, the base frequency of the imposed sine,
# The number of harmonics and a list of all the weights for separating the harmonics
sp,freq,i_har=FFT_analysis(np.array([E,i,t]).T, freq, 5, np.ones(8))
Harmonic_plots(i_har, t)
