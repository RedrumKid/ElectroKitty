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
nx=20

# Mechanism string written in full
mechanism="E(1):a+spec=b"
# Number of dissolved species, must be updated so to not brick the simulation
num_bulk_spec=2

# Constants given for adsorption and reaction on the surface
# One constant for ireversible, two for reversible
# units are as the apply, some examples are given

k0=10**2
kin_const=[
    [0.5,k0,0.],
    # [0.5,10**2,-0.1],
    # [10,1],
    # [0.5,k0,-0.1]# /, m/s, V
    # [0.5,10**4,0],                     
    # [10**7,10**-15]
    ]

# The diffusion constants for the dissolved species
# This can be any list containing the constants, must be of lenghth n_spec
# D=10**-9*np.ones(num_bulk_spec)
D=[10**-9,10**-9,10**-9]

viscosity=10**-5 #m^2/s
rot_freq=0 #Hz

# Constants describing the cell: Temperature, Ru, Cdl, A
cell_c=[T, 0, 0.*10**-2, 1*10**-4] # K, Ohm , F/m^2, m^2

# Constants describing the x direction: dx/xmax, nx
si=[0.1/36, nx, viscosity, rot_freq]

# Initial condition: first list initital surface occupations, second a list of initial condition functions, reccomend use of the uniform function 
# uniform function takes in the value of the initial concentration and give the correct way for the simulator to understand the given list

c_spec=0.01
spec_info=[
    [], #[mol/m^2]
    [uniform(c_spec/100),uniform(c_spec), uniform(0)] #[mol/m^3]
    ]

isotherm=[]

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
# time1=time.time()
# Run the simulation
E_corr,i,t=simulator_Main_loop(mechanism, [kin_const, cell_c, D, isotherm], si, t, spec_info, E) #V, A, sA
# norm=F*10**-4*F/R/T*v*10**-5
# plt.plot(t,i)
# print(max(i))
# spec, index, types=Parse_mechanism(mechanism)
# print(time.time()-time1)
# A simple plot of the calculated current
# plt.figure("Example plot", figsize=(9,5))
# plt.rc("font", size=14)
# plt.title("Example CV")
plt.plot(E,i)
# plt.ylabel("i [A]")
# plt.xlabel("E [V]")
# # plt.hlines(-0.62*F*cell_c[-1]*(D[0])**(2/3)*np.sqrt(2*np.pi*rot_freq)*(viscosity)**(-1/6),Ef,Ei,  linewidth=5, color="red")
# # plt.hlines(-0.4463,Ef,Ei, linewidth=10, color="red")
# # plt.hlines(0.25,Ef,Ei, linewidth=5, color="red")
# # plt.text(0.05,-0.3, "Randels-Sevick equation", fontsize=14)
# plt.show()

# ax.imshow(your_image, aspect='auto')
# fig.savefig(fname, dpi)

# il=-0.62*F*cell_c[-1]*(D[0])**(2/3)*np.sqrt(2*np.pi*rot_freq)*(viscosity)**(-1/6)

# ia=-(np.exp(-0.5*F/R/T*E)-np.exp(0.5*F/R/T*E))/(-np.exp(-0.5*F/R/T*E)/il+np.exp(0.5*F/R/T*E)/(-il)+1/k0/F/10**-4)

# plt.plot(E, ia)

# print("Show accuracy of simulation [%]: ")
# print((-0.62*F*cell_c[-1]*(D[0])**(2/3)*np.sqrt(2*np.pi*rot_freq)*(viscosity)**(-1/6)-min(i))/(-0.62*F*cell_c[-1]*(D[0])**(2/3)*np.sqrt(2*np.pi*rot_freq)*(viscosity)**(-1/6)))

# Functions for compting harmonics for ACV, advanced features, requieres more time points
# Uncomment for use

# sp,freq,i_har=FFT_analysis(np.array([E,i,t]).T, freq, 5, np.ones(8))
# Harmonic_plots(i_har, t)

# vs=[0.001, 0.01]

# for v in vs:
#     F=96485 #As/mol
#     R=8.314 #J/mol/K
#     T=298 #K

#     # Number of points in the spatial direction
#     nx=120

#     # Mechanism string written in full
#     mechanism="E:a=b"
#     # Number of dissolved species, must be updated so to not brick the simulation
#     num_bulk_spec=2

#     # Constants given for adsorption and reaction on the surface
#     # One constant for ireversible, two for reversible
#     # units are as the apply, some examples are given
#     kin_const=[
#         [0.5,10**-6,0.1],# /, m/s, V
#         # [0.5,10**4,0],
#         # [10**7,10**-15]
#         ]

#     # The diffusion constants for the dissolved species
#     # This can be any list containing the constants, must be of lenghth n_spec
#     # D=10**-9*np.ones(num_bulk_spec)
#     D=np.array([10**-9,1*10**-9])

#     viscosity=10**-5 #m^2/s
#     rot_freq=60 #Hz

#     # Constants describing the cell: Temperature, Ru, Cdl, A
#     cell_c=[T, 0, 0*10**-4, 10**-4] # K, Ohm , F/m^2, m^2

#     # Constants describing the x direction: dx/xmax, nx
#     si=[0.01/36, nx, viscosity, rot_freq]


#     # Initial condition: first list initital surface occupations, second a list of initial condition functions, reccomend use of the uniform function 
#     # uniform function takes in the value of the initial concentration and give the correct way for the simulator to understand the given list
#     spec_info=[
#         [], #[mol/m^2]
#         [uniform(1),uniform(0)] #[mol/m^3]
#         ]

#     # Constants used in describing the potential program
#     # In this case I use the CV programm, other are to follow, a custom one can be made
#     # by generating a the time domain and the electrical potential as numpy arrays

#     Ei=0.5  #Initial potential, also the potential to wich the CV returns [V]
#     Ef=-0.5 #Final potential, we cycle to this potential [V]
#     freq=9  #In case of ACV this is the frequency of a sine wave [Hz]
#     amp=0.  #In case of ACV this is the amplitude of a sine wave [V]
#     nt=1000 #number of time points

#     # Generate the potential program
#     E,t=V_potencial(Ei,Ef,v,amp,freq,nt,F/R/T) # V, s

#     # time1=time.time()
#     # Run the simulation
#     E,i,t=simulator_Main_loop(mechanism, [kin_const, cell_c, D], si, t, spec_info, E) #V, A, sA
#     # print(max(i))
#     # spec, index, types=Parse_mechanism(mechanism)
#     # print(time.time()-time1)
#     # A simple plot of the calculated current
#     plt.figure("Example plot", figsize=(9,5))
#     plt.title("Example CV")
#     plt.plot(E,i)
#     plt.ylabel("i [A]")
#     plt.xlabel("E [V]")




