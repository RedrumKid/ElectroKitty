# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 16:17:42 2024

@author: ozbejv
"""

############ A compact script for running a simulation
############ If you do not understand, check the tutorials

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
mechanism="E(2):a=b \n C: b=c \n C: b+*-b* \n E(1):b*=d*"

kin_const=[
    [0.5,10**2,0], # /, m/s, V
    [10,1], # /s
    [10],
    [0.5,100,-0.15]
    ]

D=3*[10**-9]

viscosity=10**-5 #m^2/s
rot_freq=0 #Hz

# Constants describing the cell: Temperature, Ru, Cdl, A
cell_c=[T, 0, 1*10**-4, 10**-4] # K, Ohm , F/m^2, m^2

# Constants describing the x direction: dx/xmax, nx
si=[0.01/36, nx, viscosity, rot_freq]

spec_info=[
    [10**-5,0,0], #[mol/m^2]
    [uniform(1),uniform(0),uniform(0)] #[mol/m^3]
    ]

isotherm=[-0.1,-0.1,-0.1]

Ei=0.3  #Initial potential, also the potential to wich the CV returns [V]
Ef=-0.3 #Final potential, we cycle to this potential [V]
v=0.1   #scan speed [V/s]
freq=9  #In case of ACV this is the frequency of a sine wave [Hz]
amp=0.0  #In case of ACV this is the amplitude of a sine wave [V]
nt=500 #number of time points

# Generate the potential program
E,t=V_potencial(Ei,Ef,v,amp,freq,nt,F/R/T) # V, s

E_corr,i,t=simulator_Main_loop(mechanism, [kin_const, cell_c, D, isotherm], si, t, spec_info, E) #V, A, sA

# A simple plot of the calculated current
plt.figure("Example plot", figsize=(9,5))
plt.title("Example CV")
plt.plot(E,i)
plt.ylabel("i [A]")
plt.xlabel("E [V]")

plt.show()

# sp,freq,i_har=FFT_analysis(np.array([E,i,t]).T, freq, 5, np.ones(8))
# Harmonic_plots(i_har, t)

