# -*- coding: utf-8 -*-
"""
Created on Wed Mar  5 14:35:00 2025

@author: ozbejv
"""

import numpy as np
from cpp_ekitty_simulator import cpp_ekitty_simulator

class gaussian_distribution:
    
    """
    A Python class for the Gaussian distribution made to work with ElectroKitty
    """
    
    def __init__(self, mu = 0, sigma = 0.1, N = 15):
        """
        Initializing by passing the mean and standard deviation
        """
        self.mu = mu
        self.sigma = sigma
        self.N = 15
        self.term = 1/np.sqrt(2*np.pi)
    
    def __call__(self, x):
        """
        The call to the class returning the probability at value x
        """
        return 1/self.sigma*self.term*np.exp(-(x-self.mu)**2/2/self.sigma**2)
    
    def update_dist_params(self, mu, sigma):
        """
        function to update mean and std 
        """
        self.mu = mu
        self.sigma = sigma
        
    def update_N(self, N):
        """
        function to update the number of points in the integration spread
        """
        self.N = N

    def create_spread(self):
        """
        function to return the area over which to integrate the parameters
        """
        return np.linspace(-3.5*self.sigma + self.mu, 3.5*self.sigma + self.mu, self.N+1)
    
    def return_mean(self):
        """
        function to return the mean
        """
        return self.mu
    
    def return_sigma(self):
        """
        function to return the sigma
        """
        return self.sigma
