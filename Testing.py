# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 09:37:03 2024

@author: ozbejv
"""
from electrokitty import ElectroKitty

F=96485 #As/mol
R=8.314 #J/mol/K
T=293 #K

# Number of points in the spatial direction
nx=20

# Mechanism string written in full
mechanism="E(1): a=b"

kin_const=[
    [0.5, 10, 0],
 ]

# The diffusion constants for the dissolved species
# This can be any list containing the constants, must be of lenghth n_spec
D=2*[10**-8]

viscosity=10**-5 #m^2/s
rot_freq=0 #Hz

# Constants describing the cell: Temperature, Ru, Cdl, A
cell_c=[T, 0, 0*1*10**-4, 10**-4] # K, Ohm , F/m^2, m^2

# Constants describing the x direction: dx/xmax, nx
si=[0.0001/36, nx, viscosity, rot_freq]

# Initial condition: first list initital surface occupations, second a list of initial condition functions, reccomend use of the uniform function 
# uniform function takes in the value of the initial concentration and give the correct way for the simulator to understand the given list
spec_info=[[1,0], [1,0]]

iso=[]
spe=[[],[1,0]]
# Constants used in describing the potential program
# In this case I use the CV programm, other are to follow, a custom one can be made
# by generating a the time domain and the electrical potential as numpy arrays

Ei=0.5 #Initial potential, also the potential to wich the CV returns [V]
Ef=-0.5 #Final potential, we cycle to this potential [V]
v=0.05   #scan speed [V/s]
freq=9 #In case of ACV this is the frequency of a sine wave [Hz]
amp=0.  #In case of ACV this is the amplitude of a sine wave [V]
nt=1000 #number of time points

problem=ElectroKitty(mechanism)
problem.V_potencial(Ei, Ef, v, amp, freq, nt)
problem.create_simulation(kin_const, cell_c, D, iso, si, spec_info)

problem.simulate()
problem.Plot_simulation()