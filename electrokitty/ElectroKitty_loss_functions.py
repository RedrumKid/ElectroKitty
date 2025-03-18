# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 12:37:58 2024

@author: ozbejv
"""

import numpy as np
import sys
import scipy.signal as scisi
from .ElectroKitty_parameter_distributions import gaussian_distribution
from cpp_ekitty_simulator import cpp_ekitty_simulator

class electrokitty_loss():
    """
    Class that calculates loss function value, given the data and the mechanism
    """
    def __init__(self, kin, species_information, cell_const, isotherm, I_data,
                 fit_kin = True,
                 fit_Cdl=False, fit_Ru=False, fit_gamamax=False,
                 fit_A=False, fit_iso=False, N_disp = 15):
        
        self.kin=kin
        self.species_information=species_information
        self.cell_const=cell_const
        self.isotherm=isotherm
        self.ysim=None
        self.I_data=I_data
        self.N_hars=5
        self.base_freq=0
        self.w=0.1*self.base_freq*np.ones(self.N_hars+1)
        self.I_har_exp=None
        self.t=None

        self.Diffusion_const = None
        self.spectators = None
        self.Spatial_info = None
        self.mechanism_list = None
        self.t = None
        self.E_generated = None
        
        self.guess, self.tells, self.gammaposition = self.create_parameter_guess(kin, species_information, 
                                                                                 cell_const, isotherm, fit_kin = fit_kin,
                                                                                 fit_Cdl=fit_Cdl, fit_Ru=fit_Ru, 
                                                                                 fit_gamamax=fit_gamamax,
                                                                                 fit_A=fit_A, fit_iso=fit_iso, N_disp=N_disp)
        
    def give_guess(self):
        """
        Gives the list of parameters that it is trying to optimise
        """
        return self.guess
    
    def set_constants(self, Diffusion_const, spectators, Spatial_info, mechanism_list, t, E):
        self.Diffusion_const = Diffusion_const
        self.spectators = spectators
        self.Spatial_info = Spatial_info
        self.mechanism_list = mechanism_list
        self.t = t
        self.E_generated = E
        
    def give_tells_gp(self):
        """
        Gives the tells list and gammaposition
        
        both are parameters that dictate the simulator how to use the guess
        """
        return self.tells, self.gammaposition
    
    def update_ysim(self, ysim):
        """
        Function that updates the function that is used when generating the simulated response
        """
        self.ysim=ysim
        
    def create_ACV_problem(self, freq, N_hars, I_har_exp, t,w):
        """
        Creating the minimisation problem for ACV
        
        needs base frequency, number of harmonics, a list containing harmonics from experimental data,
        time and a list to separate simulated harmonics
        """
        self.base_freq=freq
        self.N_hars=N_hars
        self.I_har_exp=I_har_exp
        self.t=t
        self.w=w
    
    def RMSE(self, guess):
        """
        Calculates and returns the RRMSE value given the guess for the simulator 
        """
        i_sim = self.ysim(guess)

        return np.sqrt(np.sum((i_sim-self.I_data)**2)/np.sum(self.I_data**2)/len(self.I_data))
    
    def RMSE_har(self, guess):
        """
        Calculates the harmonic average RRMSE given the guess
        """
        i_sim = self.ysim(guess)
        s, f, i_har_sim=self.FFT_analysis(self.base_freq, self.N_hars, self.w, i_sim, self.t)
        
        L=0
        for i in range(len(i_har_sim)):
            L+=1/(self.N_hars+1)*np.sqrt(np.sum((i_har_sim[i]-self.I_har_exp[i])**2)/np.sum(self.I_har_exp[i]**2)/len(self.I_har_exp[i]))
        return L
    
    def FFT_analysis(self, f,N,w, current, t):
        """
        Same as in the base class
        """
        
        def rectangular(f,w0,w):
            return np.where(abs(f-w0)<=w,1,0)
        
        I_harmonics=[]
        dt=np.average(np.diff(t))
        freq=np.fft.fftfreq(t.shape[-1],d=dt)
        sp=np.fft.fft(current)
        for i in range(N+1):
        #     #kopiram FFT
            if i==0:
                filter_sp=sp.copy()
                window=rectangular(freq,i*f,w[i])
                filter_sp=window*filter_sp
                Inew=np.fft.ifft(filter_sp).real
                I_harmonics.append(Inew)
            else:
                filter_sp=sp.copy()
                window=rectangular(freq,i*f,w[i])+rectangular(freq,-i*f,w[i])
                filter_sp=window*filter_sp
                Inew=np.real(np.fft.ifft(filter_sp))
                Inew=np.fft.ifft(filter_sp).real
                anal_signal=np.abs(scisi.hilbert(Inew))
                I_harmonics.append(anal_signal)
        return sp,freq,I_harmonics
    
    def check_type(self, param):
        """
        A simple function to check if the parameter is either a float or an intiger

        Returns:
            either False or True depending on the parameter
        """
        if type(param) is float or type(param) is int:
            return True
        else:
            return False
        
    def create_parameter_guess(self, kin, species_information, cell_const, isotherm, fit_kin=True,
                fit_Cdl=False, fit_Ru=False, fit_gamamax=False,
                fit_A=False, fit_iso=False, N_disp = 15):
        """
        Function that generates the tells list and gammamax
        
        Parameters:
            - kin: kinetic constants to be fitted
            - species_information: initial conditions
            - cell_const: cell constants
            - isotherm: isotherm constants
            
            - fit_: used to create tells
        
        tells: list containing information that tells the simulator to remap a certain parameter. 
        Contains information on which parameter it is and how to remap it. 

        tells is formed from lists, each list is made so that the parameter can be identified and placed correctly into the list it occupies.
        The first element of the list indicates the parameter kind:
            - 0: kinetic
            - 1: initial condition
            - 2: cell constant
            - 3: isotherm
        
        The second and third elements indicate the indexes to choose from in the guess. 

        In the case of cell constants the fourth element is the index in the cell_const list, 
        with the final element being a boolean, whether to fit a distribution

        In the case of initial condition the fourth element is the index in the initial condition for the adsorbed species,
        with the fith element being a boolean to fit a distribution

        In the case of isotherm the fourth element is the index of the isotherm constant, with the fifth element beeing the
        boolean to determine whether or not to fit a distribution

        In the case of kinetic constants the fourth element is the index of the list in the kinetic array, which must be
        replaced by a new list from the guess. The fifth element and onward are booleans that say which of the kinetic parameters in 
        a list is dispersed.

        gammamax: list of the lenghth of the adsorbed species initial condition and number of points to integrate
        
        """
        tells=[]
        initial_guess=[]
        
        gamma_position = [0, N_disp]

        n=0
        if fit_kin:
            for ind in range(len(kin)):
                scratch_list = [0, n, n+len(kin[ind]), ind]
                kin_guess = []
                for param in kin[ind]:
                    if self.check_type(param):
                        scratch_list.append(False)
                        kin_guess.append(param)

                    else:
                        scratch_list.append(True)
                        scratch_list[2] += 1
                        kin_guess.append(param[0].return_mean())
                        kin_guess.append(param[0].return_sigma())

                n = scratch_list[2]
                tells.append(scratch_list)
                
                initial_guess+=kin_guess

        if fit_Ru:
            if self.check_type(cell_const[1]):
                tells.append([2, n, n+1, 1, False])
                initial_guess.append(cell_const[1])
                n+=1
            else:
                tells.append([2, n, n+2, 1, True])
                initial_guess.append(cell_const[1][0].return_mean())
                initial_guess.append(cell_const[1][0].return_sigma())
                n+=2

        if fit_Cdl:
            if self.check_type(cell_const[2]):
                tells.append([2, n, n+1, 2, False])
                initial_guess.append(cell_const[2])
                n+=1
            else:
                tells.append([2, n, n+2, 2, True])
                initial_guess.append(cell_const[2][0].return_mean())
                initial_guess.append(cell_const[2][0].return_sigma())
                n+=2
            
        if fit_A:
            if self.check_type(cell_const[3]):
                tells.append([2, n, n+1, 3, False])
                initial_guess.append(cell_const[3])
                n+=1
            else:
                tells.append([2, n, n+2, 3, True])
                initial_guess.append(cell_const[3][0].return_mean())
                initial_guess.append(cell_const[3][0].return_sigma())
                n+=2
        
        if fit_gamamax:
            gama_guess = []
            gamma_position[0] = len(species_information[0])

            for ind in range(len(species_information[0])):
                if self.check_type(species_information[0][ind]) and species_information[0][ind]>0:
                    scratch_list = [1, n, n+1, ind, False]
                    gama_guess.append(species_information[0][ind])
                    tells.append(scratch_list)
                    n += 1
                elif self.check_type(species_information[0][ind]) == False:
                    scratch_list = [1, n, n+2, ind, True]
                    gama_guess.append(species_information[0][ind][0].return_mean())
                    gama_guess.append(species_information[0][ind][0].return_sigma())
                    tells.append(scratch_list)
                    n += 2
            
            initial_guess+=gama_guess
        
        if fit_iso:
            
            iso_guess = []

            for ind in range(len(isotherm)):
                if self.check_type(isotherm[ind]):
                    scratch_list = [3, n, n+1, ind, False]
                    iso_guess.append(isotherm[ind])
                    tells.append(scratch_list)
                    n += 1
                else:
                    scratch_list = [3, n, n+2, ind, True]
                    iso_guess.append(isotherm[ind][0].return_mean())
                    iso_guess.append(isotherm[ind][0].return_sigma())
                    tells.append(scratch_list)
                    n += 2
            
            initial_guess+=iso_guess
        
        return np.array(initial_guess), tells, gamma_position
    
    def unpack_fit_params(self, guess, tells, gamma_position, kin, species_information, cell_const, isotherm):
        """
        Function takes the guess, tells and gammma_position to reconstruct the lists for the simulator

        Parameters:
            - guess: the list of parameters that are beeing fitted
            - tells: the list with information on how to remap the guess into correct lists
            - gamma_position: a list that has helpful information for remaping
            - kin: the list containing kinetic parameters
            - species_information: the list with the initial condition
            - cell_const: the list with the cell constants
            - isotherm: the list with isotherm constants
        
        Returns:
            - kin: the remapped kinetic constants
            - cell_const: the remapped cell constants
            - species_information: the remapped initial condition
            - isotherm: the remapped isotherm constants
        """
        guess=guess.tolist()
        kinetics=list(kin)
        cell_params=list(cell_const)
        spec_info=list(species_information)
        iso = list(isotherm)
        
        for info in tells:
            if info[0] == 0:
                temp = guess[info[1]:info[2]]
                kin_list = []
                count = 0
                for ind in range(len(info[4:])):
                    if info[4+ind]:
                        if len(info[4:]) == 3 and ind == 1:
                            kin_list.append([gaussian_distribution(temp[count], temp[count+1]), gamma_position[1],
                                            -3.5*temp[count+1]+temp[count], 3.5*temp[count+1]+temp[count], "log"])
                        else:
                            kin_list.append([gaussian_distribution(temp[count], temp[count+1]), gamma_position[1],
                                            -3.5*temp[count+1]+temp[count], 3.5*temp[count+1]+temp[count], "lin"])
                        count += 2
                    else:
                        kin_list.append(temp[count])
                        count += 1
                kinetics[info[3]] = kin_list
            
            elif info[0] == 1:
                if info[2]-info[1] == 1 and info[-1] == False:
                    spec_info[0][info[3]] = guess[info[1]]
                else:
                    spec_info[0][info[3]] = [gaussian_distribution(guess[info[1]], guess[info[1]+1]), gamma_position[1],
                                            -3.5*guess[info[2]-1]+guess[info[1]],
                                            3.5*guess[info[2]-1]+guess[info[1]], "lin"]
            elif info[0] == 2:
                if info[2]-info[1] == 1 and info[-1] == False:
                    cell_params[info[3]] = guess[info[1]]
                else:
                    cell_params[info[3]] = [gaussian_distribution(guess[info[1]], guess[info[1]+1]), gamma_position[1],
                                            -3.5*guess[info[2]-1]+guess[info[1]],
                                            3.5*guess[info[2]-1]+guess[info[1]], "lin"]

            elif info[0] == 3:
                if info[2]-info[1] == 1 and info[-1] == False:
                    iso[info[3]] = guess[info[1]]
                else:
                    iso[info[3]] = [gaussian_distribution(guess[info[1]], guess[info[1]+1]), gamma_position[1],
                                            -3.5*guess[info[2]-1]+guess[info[1]],
                                            3.5*guess[info[2]-1]+guess[info[1]], "lin"]


        return kinetics, cell_params, spec_info, iso

    
    def create_lower_upper_bounds(self, guess, tells, potential):
        """
        Function creates based on the guess (arameters to be fitted) the lower and upper bounds,
        for either CMA-ES or MCMC. Both algorithms follow the same rules
        
        The default is:
            - alpha: 0,1
            - k0: 0, 100*k0 or 1000
            - E0: -0.5+min(E), 0.5+max(E)
            - kf, kb: 0, 100*k
            - Cdl: 0, 100*Cdl
            - Ru: 0, 100*Ru
            - A: 0, 100*A
            - gammamax: 0, 100*gammamax
            - isotherm: -25, 10

        Parameters:
            - guess: list of parameters that will be fitted
            - tells: the list containing information to remap the guess into correct simulation lists
            - potential: the potential program used in the simulation
        
        Returns:
            - the lower bounds for constrained optimisation or MCMC sampling
            - the upper bounds for constrained optimisation or MCMC sampling
        """
        guess=guess.tolist()
        lower_bound=[]
        upper_bound=[]
        
        pot_min, pot_max = min(potential)-0.5, max(potential)+0.5
        for info in tells:
            if info[0] == 0:
                temp = guess[info[1]:info[2]]
                count = 0
                if len(info[4:]) == 3:
                    for ind in range(len(info[4:])):
                        if ind == 0:
                            if info[4+ind]:
                                lower_bound.append(0)
                                lower_bound.append(10**-4)
                                upper_bound.append(1)
                                if temp[count+1] == 0:
                                    upper_bound.append(1)
                                else:
                                    upper_bound.append(100*temp[count+1])
                                count += 2
                            else:
                                lower_bound.append(0)
                                upper_bound.append(1)
                                count += 1

                        elif ind == 1:
                            if info[4+ind]:
                                lower_bound.append(0)
                                lower_bound.append(10**-4)
                                if temp[count] == 0:
                                    upper_bound.append(1000)
                                else:
                                    upper_bound.append(100*temp[count])

                                if temp[count+1] == 0:
                                    upper_bound.append(1)
                                else:
                                    upper_bound.append(100*temp[count+1])
                                count += 2
                            else:
                                lower_bound.append(0)
                                if temp[count] == 0:
                                    upper_bound.append(1000)
                                else:
                                    upper_bound.append(100*temp[count])
                                count += 1

                        elif ind == 2:
                            if info[4+ind]:
                                lower_bound.append(pot_min)
                                lower_bound.append(10**-4)
                                upper_bound.append(pot_max)
                                upper_bound.append(100*temp[count+1])
                                count += 2
                            else:
                                lower_bound.append(pot_min)
                                upper_bound.append(pot_max)
                                count += 1
                elif len(info[4:]) == 2 or len(info[4:]) == 1:
                    for ind in range(len(info[4:])):
                        if info[4+ind]:
                                lower_bound.append(0)
                                lower_bound.append(10**-4)
                                if temp[count] == 0:
                                    upper_bound.append(1000)
                                else:
                                    upper_bound.append(100*temp[count])

                                if temp[count+1] == 0:
                                    upper_bound.append(1)
                                else:
                                    upper_bound.append(100*temp[count+1])
                                count += 2
                        else:
                            lower_bound.append(0)
                            if temp[count] == 0:
                                    upper_bound.append(1000)
                            else:
                                upper_bound.append(100*temp[count])
                            count += 1

            if info[0] != 0 and info[0] != 3:
                if info[-1]:
                    lower_bound.append(0)
                    lower_bound.append(10**-4)
                    if guess[info[1]] == 0:
                        if info[0] == 1:
                            upper_bound.append(1)
                        elif info[0] == 2 and info[3] == 1:
                            upper_bound.append(100)
                        elif info[0] == 2 and info[3] == 2:
                            upper_bound.append(100)
                    else:
                        upper_bound.append(100*guess[info[1]])
                    if guess[info[1]+1] == 0:
                        upper_bound.append(1)
                    else:
                        upper_bound.append(100*guess[info[1]+1])
                else:
                    lower_bound.append(0)
                    if guess[info[1]] == 0:
                        if info[0] == 1:
                            upper_bound.append(1)
                        elif info[0] == 2 and info[3] == 1:
                            upper_bound.append(100)
                        elif info[0] == 2 and info[3] == 2:
                            upper_bound.append(100)
                    else:
                        upper_bound.append(100*guess[info[1]])

            if info[0] == 3:
                if info[-1]:
                    lower_bound.append(-25)
                    lower_bound.append(10**-4)
                    
                    upper_bound.append(10)

                    if guess[info[1]+1] == 0:
                        upper_bound.append(1)
                    else:
                        upper_bound.append(100*guess[info[1]+1])

                else:
                    lower_bound.append(-25)
                    upper_bound.append(10)

        lower_bound.append(0)
        upper_bound.append(1)

        return lower_bound, upper_bound 
    
    def create_axis_labels_old(self, tells, a_spec):
        """
        Function takes the tells and species names to return parameter labels for the chain

        !!! This is the old version and is used solely because of previous saved files from electrokittys

        Parameters:
            - tells: the list containing remaps from the guess into correct simulation lists
            - a_spec: a list containing the names of all adsorbed species
        
        Returns:
            - the list of labels for the parameters sampled during MCMC
        """
        
        labels = []
        index1 = 0
        for i in range(tells[0]):
                
            index2=tells[i+1]
            dist = abs(index2-index1)
            if dist == 3:
                labels.append(r"$\alpha"+"_{"+str(i+1)+"}$")
                labels.append(r"$k_{0,"+str(i+1)+"}$")
                labels.append(r"$E^{0'}_{"+str(i+1)+"}$")
            elif dist == 2:
                labels.append(r"$k_{f"+str(i+1)+"}$")
                labels.append(r"$k_{b"+str(i+1)+"}$")
            elif dist == 1:
                labels.append(r"$k_{"+str(i+1)+"}$")
            index1=index2

        if tells[tells[0]+1] != 0:
            labels.append(r"$R_{u}$") #Ru
        else:
            pass
        
        if tells[tells[0]+2] != 0:
            labels.append(r"$C_{dl}$") #Cdl
        else:
            pass
        
        if tells[tells[0]+3] != 0:
            labels.append(r"$A$") #A
        else:
            pass
        
        if tells[tells[0]+4] != 0:
            labels.append(r"$\Gamma_{max}$") #gammamax
        else:
           pass
           
        if tells[tells[0]+5]!=0:
            for spec_name in a_spec: #isotherm
                labels.append(r"$g_{"+spec_name+"}$")
        else:
            pass
        labels.append(r"$\sigma$")
        return labels
    
    def create_axis_labels(self, tells, a_spec):
        """
        Function takes the tells and species names to return parameter labels for the chain

        Parameters:
            - tells: the list containing remaps from the guess into correct simulation lists
            - a_spec: a list containing the names of all adsorbed species
        
        Returns:
            - the list of labels for the parameters sampled during MCMC
        """
        labels = []
        for index in range(len(tells)):
            info = tells[index]
            if info[0] == 0:
                count = 0
                if len(info[4:]) == 3:
                    for ind in range(len(info[4:])):
                        if ind == 0:
                            if info[4+ind]:
                                labels.append(r"$\alpha"+"_{"+str(index+1)+"}$")
                                labels.append(r"\sigma_{\alpha}")
                                count += 2
                            else:
                                labels.append(r"$\alpha"+"_{"+str(index+1)+"}$")
                                count += 1

                        elif ind == 1:
                            if info[4+ind]:
                                labels.append(r"$k_{0,"+str(index+1)+"}$")
                                labels.append(r"$\sigma_{k0}$")
                                count += 2
                            else:
                                labels.append(r"$k_{0,"+str(index+1)+"}$")
                                count += 1

                        elif ind == 2:
                            if info[4+ind]:
                                labels.append(r"$E^{0'}_{"+str(index+1)+"}$")
                                labels.append(r"$\sigma_{E0}$")
                                count += 2
                            else:
                                labels.append(r"$E^{0'}_{"+str(index+1)+"}$")
                                count += 1

                elif len(info[4:]) == 2:
                    for ind in range(len(info[4:])):
                        if info[4+ind]:
                            labels.append(r"$k_{"+str(index+1)+"}$")
                            labels.append(r"$\sigma_{k}")
                            count += 2
                        else:
                            labels.append(r"$k_{"+str(index+1)+"}$")
                            count += 1
                
                elif len(info[4:]) == 1:
                    for ind in range(len(info[4:])):
                        if info[4+ind]:
                            labels.append(r"$k_{"+str(index+1)+"}$")
                            labels.append(r"$\sigma_{k}")
                            count += 2
                        else:
                            labels.append(r"$k_{"+str(index+1)+"}$")
                            count += 1

            if info[0] == 1:
                if info[-1]:
                    labels.append(r"$\Gamma_{"+a_spec[info[3]]+"}$")
                    labels.append(r"$\sigma_{\Gamma}$")
                else:
                    labels.append(r"$\Gamma_{"+a_spec[info[3]]+"}$")
            
            if info[0] == 2:
                if info[3] == 1:
                    if info[-1]:
                        labels.append(r"$R_{u}$")
                        labels.append(r"$\sigma_{Ru}$")
                    else:
                        labels.append(r"$R_{u}$")

                if info[3] == 2:
                    if info[-1]:
                        labels.append(r"$C_{dl}$")
                        labels.append(r"$\sigma_{Cdl}$")
                    else:
                        labels.append(r"$C_{dl}$")

                if info[3] == 3:
                    if info[-1]:
                        labels.append(r"$A$")
                        labels.append(r"$\sigma_{A}$")
                    else:
                        labels.append(r"$A$")

            if info[0] == 3:
                if info[-1]:
                    labels.append(r"$g_{"+a_spec[info[3]]+"}$")
                    labels.append(r"$\sigma_{g}$")
                else:
                    labels.append(r"$g_{"+a_spec[info[3]]+"}$")

        labels.append(r"$\sigma_{exp}$")
        return labels
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        