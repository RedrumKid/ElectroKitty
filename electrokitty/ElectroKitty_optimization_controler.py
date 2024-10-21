# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 11:18:50 2024

@author: ozbejv
"""

import numpy as np
from .optimization_routines import electrokitty_optimization_algorithms

class electrokitty_optimization_controller:
    """
    Wrapper class that correctly executes the algorithm for parameter estimation
    """
    def __init__(self,fun, x0,
    algorithm="Nelder-Mead", tolf = 10**-11, tolx = 10**-11, lb=None, ub=None):
        
        self.tells=None
        self.gammaposition=None
        self.best_guess=None
        self.fun=fun
        self.ub=ub
        self.lb=lb
        self.tolf=tolf
        self.tolx=tolx
        self.fit_score=None
        self.algorithm=algorithm
        self.x0=x0
        
        self.optimizer=electrokitty_optimization_algorithms()
    
    def fit_parameters(self):
        
        if self.algorithm=="Nelder-Mead":
            sol, message= self.optimizer.nelder_mead(self.fun, self.x0, no_improve_thershold=self.tolf)
            self.best_guess=sol[0]
            self.fit_score=sol[1]
            print("Final message from Nelder-Mead:")
            print(message)
            print()
        
        elif self.algorithm=="CMA-ES-inbuilt":
            sol, message, succ=self.optimizer.CMA_controller(self.fun, self.x0, 
                                                    0.1*np.ones(len(self.x0)),
                                                   lb=self.lb, ub=self.ub,
                                                   tolf=self.tolf, tolx=self.tolx)
            self.best_guess=sol[0]
            self.fit_score=sol[1]
            print()
            print("Final message from CMA-ES:")
            print(message)
            print()
            if succ:
                print("CMA-ES found solution")
            else:
                print("CMA-ES failed")
        
        elif self.algorithm=="CMA-ES":
            
            sol, message = self.optimizer.CMA_ES(self.fun, self.x0, 0.1, 
                                                 lb=self.lb, ub=self.ub, tolf=self.tolf, tolx=self.tolx)
            
            self.best_guess=sol[0]
            self.fit_score=sol[1]
        
        return self.best_guess, self.fit_score
    
            