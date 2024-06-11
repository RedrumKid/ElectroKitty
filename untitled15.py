# -*- coding: utf-8 -*-
"""
Created on Tue May 21 09:49:53 2024

@author: ozbejv
"""

import numpy as np
import matplotlib.pyplot as plt

def f_i(E, k0, D):
    
    kf=k0*np.exp(-0.5*38.9*E)
    kb=k0*np.exp(0.5*38.9*E)
    
    i=1/(96500*10**-4*(kf-kb))*(1+(D**(-2/3)*kf+D**(-2/3)*kb)/(0.62*(10**-5)**(-1/6)*np.sqrt(2*np.pi*60)))
    return 1/i

E=np.linspace(-0.5,0.5,100)

i=f_i(E,100000,10**-9)

plt.plot(E,i)