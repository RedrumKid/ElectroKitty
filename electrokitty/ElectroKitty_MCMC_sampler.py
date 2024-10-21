# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 14:58:33 2024

@author: ozbejv
"""

import numpy as np
import multiprocessing

class electrokitty_sampler():
    """
    Class containing the MCMC sampler used in parameter estimation
    """
    def __init__(self, n_samples, burn_in_per, num_chains, 
                 multi_processing, bounds, I_data):
        
        self.n_samples=n_samples
        self.burn_in_per=burn_in_per
        self.num_chains=num_chains
        self.multi_processing=multi_processing
        self.bounds=bounds
        self.y_sim=None
        self.I_data=I_data
    
    def give_y_sim(self, ysim):
        """
        update the function used to generate the simulation
        """
        self.y_sim=ysim
    
    def prior(self, theta, lower_bound, upper_bound):
        """
        Function that given the guess and the bounds calculates the prior,
        either 0 or 1
        """
        p=[]
        
        for i in range(len(theta)):
            
            if theta[i] >= lower_bound[i] and theta[i] <= upper_bound[i]:
                p.append(1)
            else:
                p.append(-np.inf)
        return min(p)
    
    def likelihood(self, theta):
        """
        Given the guess calculates the likelihood based on a guessian that the guess fits the data
        """
        i_sim=self.y_sim(theta)
        
        N=len(self.I_data)
            
        p=-N/2*np.log(2*np.pi*theta[-1]**2)-np.sum((self.I_data-i_sim)**2)/2/theta[-1]**2
            
        return p
    
    def propose(self, p, covar):
        """
        A function that given the means and their covariance matrix proposes a guess
        """
        
        p_new=np.random.multivariate_normal(p, covar)

        return p_new           

    def accept(self, p, p_new, lower_bound, upper_bound, acc_rate):
        """
        Calculates the likelihood of the old and new guess and either accepts it or rejects the new guesss
        returning either
        """
        theta=p
        acc=0
        
        numerator = self.likelihood(p) + (self.prior(p, lower_bound, upper_bound))
        
        numerator_new = self.likelihood(p_new) + (self.prior(p_new, lower_bound, upper_bound))
        
        if numerator_new >= numerator:
            theta=p_new
            acc_rate+=1
            acc=1

        elif np.log(np.random.rand()) < (numerator_new-numerator):
            theta=p_new
            acc_rate+=1
            acc=1
            
        return theta, acc_rate, acc
    
    def update_moms(self, means, M2, current, log_lam, i, p):
        """
        updates the parameter mean and covariance matrix based on the new guess
        
        using HB-MCMC it updates the mean and covariance based on the acceptance rate
        """
        nexti=i+1
        w=nexti**-0.5
        newme=(1-w)*means+w*current
        db, da = (current-newme), (current-newme)
        newM=(1-w)*M2+w*np.outer(db,da)
        log_lam+=w*(p-0.225)
        return newme, newM, log_lam
    
    def single_chain(self, n_samples, means, lower_bound, upper_bound, chain_num):
        """
        Calculates a single MCMC chain using HBMCMC
        """
        divisor=max([1,int(0.01*n_samples)])
        results=[]
        
        dim=len(means)
        
        cov0=1/np.sqrt(dim)*np.eye(dim)
        
        results.append(means)
        
        M2=cov0
        means=means
        acc=0
        log_lam=0
        acc_total=0
        
        for i in range(n_samples):
            x_cur=results[-1]
            
            prop_theta=self.propose(x_cur, np.exp(log_lam)*M2)
            
            theta, acc_total, acc = self.accept(x_cur, prop_theta, lower_bound, upper_bound, acc_total)
            results.append(theta)
            
            means, M2, log_lam = self.update_moms(means, M2, results[-1], log_lam, i, acc_total/(i+1))
            if i%divisor==0:
                print(f"chain num: {chain_num}   iteration: {i}   acceptance rate: {round(acc_total/(i+1),3)}")
            
        return np.array(results)
    
    def single_chain_wrap(self, args):
        """
        Call to execute a single chain
        """
        return self.single_chain(*args)
    
    def start_sampler(self, initial_guess):
        """
        Call to start samplig from the initial guess
        
        This will perform either on a single or multicore depending on what was chosen
        """
        lower_bound, upper_bound=self.bounds
        
        initial_positions=[]
        scatter=np.random.uniform(0.8, 1.2, self.num_chains)
        for index in range(self.num_chains):
            initial_positions.append(scatter[index]*initial_guess)
        
        chains=[]
        
        if self.multi_processing == False:
            print("Starting MCMC sampling")
            print()
            
            for ch_num in range(self.num_chains):
                print()
                print(f"Resolving chain num: {ch_num}")
                print()
                single_chain=self.single_chain(self.n_samples,
                                                initial_positions[ch_num],
                                                    lower_bound, upper_bound, ch_num)
                chains.append(single_chain)
            
            
        else:
            stuff=[]
            for i in range(self.num_chains):
                stuff.append([self.n_samples, initial_positions[i], 
                              lower_bound, upper_bound, i])
            
            pool=multiprocessing.Pool(self.num_chains)
            chains=pool.map(self.single_chain_wrap, stuff)
            pool.close()
            pool.join()
            
        
        chains=np.array(chains)
        return chains
    
    