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
T=293 #K

# Number of points in the spatial direction
nx=20

# Mechanism string written in full
mechanism="E(2):a=b \n C: b=c \n C: b+*-b* \n E(1):b*=d*"
# Number of dissolved species, must be updated so to not brick the simulation
num_bulk_spec=3

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
D=10**-9*np.ones(num_bulk_spec)

viscosity=10**-5 #m^2/s
rot_freq=0 #Hz

# Constants describing the cell: Temperature, Ru, Cdl, A
cell_c=[T, 30, 1*10**-4, 10**-4] # K, Ohm , F/m^2, m^2

# Constants describing the x direction: dx/xmax, nx
si=[0.1/36, nx, viscosity, rot_freq]


# Initial condition: first list initital surface occupations, second a list of initial condition functions, reccomend use of the uniform function 
# uniform function takes in the value of the initial concentration and give the correct way for the simulator to understand the given list
spec_info=[
    [10**-5,0,0], #[mol/m^2]
    [uniform(1),uniform(0),uniform(0)] #[mol/m^3]
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

# time1=time.time()
# Run the simulation
E,i,t=simulator_Main_loop(mechanism, [kin_const, cell_c, D], si, t, spec_info, E) #V, A, sA
# spec, index, types=Parse_mechanism(mechanism)
# print(time.time()-time1)
# A simple plot of the calculated current
plt.figure("Example plot", figsize=(9,5))
plt.title("Example CV")
plt.plot(E,i)
plt.ylabel("i [A]")
plt.xlabel("E [V]")
# plt.hlines(-0.62*F*cell_c[-1]*(D[0])**(2/3)*np.sqrt(2*np.pi*rot_freq)*(viscosity)**(-1/6),Ef,Ei)
plt.show()

# print("Show accuracy of simulation [%]: ")
# print((-0.62*F*cell_c[-1]*(D[0])**(2/3)*np.sqrt(2*np.pi*rot_freq)*(viscosity)**(-1/6)-min(i))/(-0.62*F*cell_c[-1]*(D[0])**(2/3)*np.sqrt(2*np.pi*rot_freq)*(viscosity)**(-1/6)))

# Functions for compting harmonics for ACV, advanced features, requieres more time points
# Uncomment for use

# sp,freq,i_har=FFT_analysis(np.array([E,i,t]).T, freq, 5, np.ones(8))
# Harmonic_plots(i_har, t)

# A way to simulate impedance 
# Must still be chacked and verified
# Uncomment to use, note that it is slow

# freqs=np.logspace(-3,4,40)
# currents=[]
# potentials=[]
# for freq in freqs:
#     # freq=10**fr
#     t=np.linspace(0,10/freq,201)
#     E=0.005*np.sin(2*np.pi*freq*t)-0.
    
#     dt=t[1]
    
#     mechanism="E:a=b"

#     Ads_c=[
#             # [10]
#             ]

#     B_c=[
#            # [10,1]
#           ]

#     ec_c=[
#           [0.5,0*10**2,0],
#             # [0.5,100,0]
#           ]

#     D=10**-9*np.ones(2)

#     cell_c=[293, 10, 1*10**-4, 10**-4]

#     si=[0.1/36,15]

#     spec_info=[
#         [],
#         [uniform(0.5), uniform(0.5)]
#         ]

#     # E,t=V_potencial(0.3,-0.3,0.1,0.1,9,10000,96485/8.314/293)

#     E,i,t=simulator_Main_loop(mechanism, [Ads_c, B_c, ec_c, cell_c, D], si, t, spec_info, E)
#     currents.append(i)
#     potentials.append(E)
    
# r=[]
# im=[]
# for i in range(len(freqs)):
#     j=np.fft.fft(currents[i])
#     e=np.fft.fft(potentials[i])
#     impedance=j/e
#     # impedance=0.005/j
#     # plt.plot(e)
#     # plt.plot(j)
#     # j=j[10]
#     # e=e[10]
#     # plt.plot(impedance[2:])
#     r.append(impedance[10].real)
#     im.append(impedance[10].imag)

# r=np.array(r)
# im=np.array(im)
# plt.plot(r,im)













