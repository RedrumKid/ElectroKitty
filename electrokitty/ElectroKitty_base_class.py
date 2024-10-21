# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 14:43:51 2024

@author: ozbejvodeb
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
from datetime import datetime
import scipy.optimize as sciop
import scipy.signal as scisi
import _pickle as cPickle
import webbrowser
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from .ElectroKitty_parser import electrokitty_parser
from .ElectroKitty_loss_functions import electrokitty_loss
from .ElectroKitty_optimization_controler import electrokitty_optimization_controller
from .ElectroKitty_MCMC_sampler import electrokitty_sampler
from cpp_ekitty_simulator import cpp_ekitty_simulator

class ElectroKitty:
    """
    ElectroKitty object for simulation and data fitting
    Includes ways of making 2 potential programs for CV/ACV and chronoamperometry
    It computes the current response for an RDE setup, with or without rotation
    You can plot the simulation, data, adsorbed species profile and concentration profile
    For data fitting available algorithms are Nelder-Mead and CMA-ES
    Also included is a MCMC routine to sample a parameter distribution
    """
    
    def __init__(self, string="E(1):a=b"):
        """
        initialise parameters
        """
        self.F=96485
        self.R=8.314
        self.t=None
        self.sp=None
        self.freq=None
        self.I_harmonics=None  
        self.I_har_data=None
        self.E_generated=None
        self.current=None
        self.I_data=None
        self.concentration_profile=None
        self.surface_profile=None
        
        self.string=string
        
        self.cell_const=None
        self.diffusion_const=None
        self.isotherm=None
        self.spectators=None
        self.spatial_info=None
        self.species_information=None
        self.kin=None
        
        self.x=None
        self.number_of_diss_spec=None
        self.number_of_surf_conf=None
        self.E_Corr=None
        self.mechanism_list=None
        
        self.fit_score=None
        self.tells=None
        self.gamaposition=None
        self.multi_core_MCMC=False
        self.chains=None
        self.mean_chain=None
        
        self.Parser=electrokitty_parser(string)
        self.simulator=cpp_ekitty_simulator()
        self.optimizer=None
        self.loss_function=None
        self.MCMC_sampler=None
    
    def create_simulation(self, kin, cell_const, 
                          Diffusion_const, isotherm,Spatial_info, 
                          Species_information, spectators=False):
        """
        Give the simulator all required parameters for simulation
        
        The parameters are: 
            kin: kinetic constants to simulate from, also used as a first guess when fitting
                EC: kinetics need a list of 2, C: needs either 2 or 1, depending on the mechanism
            
            cell_const: list of 4 constants in order of T, Ru, Cdl, A
            
            Diffusion_const: a list containing all diffusion constants for all dissolved species
            
            isotherm: a list of isotherm constants for all adsorbed species
            
            Spatial_info: a list of constants for creating the diffusion layer in order:
                dx/xmax, nx (number of points), viscosity, rotational frequency in Hz
            
            Species_information: initial condition of the simulation,
                a list of lists: surface concentrations, dissolved concentrations
            
            spectators: not fully implemented, a list of 1 or 0 to use the species or not in the reaction
        
        """
        
        self.cell_const=cell_const
        self.diffusion_const=Diffusion_const
        self.number_of_surf_conf=len(Species_information[0])
        self.number_of_diss_spec=len(Species_information[1])
        self.isotherm=isotherm
        self.spectators=spectators
        self.spatial_info=Spatial_info
        self.species_information=Species_information
        self.kin=kin
        
        spectators = [np.ones(len(Species_information[0])),np.ones(len(Species_information[1]))]
        self.spectators = spectators

        self.mechanism_list=self.Parser.Parse_mechanism()
        self.simulator.set_parameters(
                              cell_const, Diffusion_const, isotherm, spectators, Spatial_info, Species_information, kin, 
                              self.mechanism_list[0], self.mechanism_list[1], 
                              self.mechanism_list[2], self.mechanism_list[3], self.mechanism_list[4]
                              
                              )

        self.simulator.set_simulation_programm(self.t, self.E_generated)
	
    def save(self, filename):
        """
        function to save all class parameters in a file.
        
        the filename is the file's name, .ek is appended inside
        !!! The file is not human readable, use load function to update a class
        """
        
        f=open(filename+".ek", "wb")
        
        save_list = [
            self.t,
            self.sp,
            self.freq,
            self.I_harmonics,   
            self.E_generated,
            self.current,
            self.I_data,
            self.concentration_profile,
            self.surface_profile,
            self.string,
            self.cell_const,
            self.diffusion_const,
            self.isotherm,
            self.spectators,
            self.spatial_info,
            self.species_information,
            self.kin,
            self.x,
            self.number_of_diss_spec,
            self.number_of_surf_conf,
            self.E_Corr,
            self.mechanism_list,
            self.fit_score,
            self.tells,
            self.gamaposition,
            self.multi_core_MCMC,
            self.chains,
            self.mean_chain
            ]
        
        cPickle.dump(save_list, f, 2)
        f.close()
    
    def load(self, filename):
        """
        a function to load a .ek file containing ElectroKitty class parameters to a class
        the function loads and updates an ElectroKitty class instantce
        """
        f=open(filename, "rb")
        temp_dict = cPickle.load(f)
        f.close()
        self.t, self.sp, self.freq, self.I_harmonics, self.E_generated = temp_dict[:5]
        self.current, self.I_data, self.concentration_profile, self.surface_profile  = temp_dict[5:9]
        self.string ,self.cell_const, self.diffusion_const, self.isotherm, self.spectators = temp_dict[9:14]
        self.spatial_info, self.species_information, self.kin, self.x = temp_dict[14:18]
        self.number_of_diss_spec, self.number_of_surf_conf, self.E_Corr = temp_dict[18:21]
        self.mechanism_list, self.fit_score, self.tells, self.gamaposition = temp_dict[21:25]
        self.multi_core_MCMC,self.chains,self.mean_chain = temp_dict[25:]
        
        self.Parser=electrokitty_parser(self.string)
        self.mechanism_list=self.Parser.Parse_mechanism()
        self.simulator.set_parameters(
                              self.cell_const, self.diffusion_const, self.isotherm, self.spectators, 
                              self.spatial_info, self.species_information, self.kin, 
                              self.mechanism_list[0], self.mechanism_list[1], 
                              self.mechanism_list[2], self.mechanism_list[3], self.mechanism_list[4]
                              )

        self.simulator.set_simulation_programm(self.t, self.E_generated)
    
    def set_data(self, E_data, i_data, t_data):
        """
        a function for importing data. This function updates E_generated the potential signal used for simulation
        important when fitting, as well as when trying to simulate current to custom built signals
        """
        self.I_data=i_data
        self.E_generated=E_data
        self.t=t_data 
    
    def load_data_from_txt(self, path, skip=0):
        """
        Function for loading data from .txt file

        Must be formatted, E, i, t
        Delimiter is a tab

        Parameters:
            pats... the path to the file
            skip... skips this many rows in the data
        """
        data = np.loadtxt(path, skiprows=skip)
        self.I_data = data[:,1]
        self.E_generated = data[:,0]
        self.t = data[:,2]
    
    def IR_correct_data(self, Ru=-1, correct_data = False, correct_simulation = False):
        """
        Function for correcting the Ru in data or simulation

        Parameters:
            Ru... the ammount of resistance to be corrected [Ohm]
            must choose either correct_data or correct_simulation
        """
        if Ru < 0:
            Ru = self.cell_const[1]
        if correct_simulation:
            self.E_corrected = self.E_generated - Ru*self.current
        elif correct_data:
            self.E_corrected = self.E_generated - Ru*self.I_data
        else:
            print("Nothing to correct: choose either data or simulation")
    
    def simulate(self, eqilib=False):
        """
        function to call the simulator.
        
        Will update the current, E_Corr, surface_profile, concentration_profile
        """
        
        self.current = self.simulator.simulate()
        self.E_Corr = self.simulator.give_E_corr()
        self.surface_profile = self.simulator.give_surf_profile()
        self.concentration_profile = self.simulator.give_concentration_profile()
        # self.E_Corr, self.current, self.surface_profile, self.concentration_profile = self.simulator.simulate(eqilib=eqilib)
        
    ####################### Functions for generating potential programs
    
    def V_potencial(self, Ei,Ef,v,amplitude,frequency,nt):
        """
        Function to generate a potential signal for a CV/ACV simulation
        
        Parameters:
            Ei... initial potential [V]
            Ef... potential when the signal turns [V]
            v... scan rate [V/s]
            amplitude... sinusoidal amplitude [V]
            frequency... sine frequency [Hz]
            nt... number of time points
        
        !!! The function will override the data potential signal
        """
        # calculate the potentail input signal for a CV or ACV
        En=(Ei+Ef)/2
        tp=2*(Ei-Ef)/v
        a=Ei-En
        ts=abs(2*(-Ei+Ef)/v)
        self.t=np.linspace(0,ts,nt)
        self.E_generated=a-2*a/np.pi*np.arccos(np.cos(2*np.pi/tp*self.t))+En+amplitude*np.sin(2*np.pi*frequency*self.t)
        return self.E_generated,self.t
    
    def C_potential(self, E, tmax, nt):
        """
        Function to generate a constant potential for chronoamperometry
        
        Parameters:
            E... holding potential [V]
            tmax... maximum time [s]
            nt... number of time points
        """
        self.E_generated=E*np.ones(nt)
        self.t=np.linspace(0,tmax,nt)
        return self.E_generated, self.t
    ##################### Fitting to data
    
    
    def fit_to_data(self, fit_Cdl=False, fit_Ru=False, fit_gamamax=False,
                    fit_A=False, fit_iso=False, eqilibration=False, algorithm="Nelder-Mead",
                    tolf=10**-11, tolx=10**-11):
        
        """
        Function which tries to fit kinetic parameters and others to current from data
        
        Parameters:
            fit_Cdl... True to fit the double layer capacitance
            fit_Ru... True to fit uncompensated resistance
            fit_gamamax... True to fit the highest surface concentration in the species_information list
            fit_A... True to fit the electrode surface
            fit_iso... True to fit the whole list of isotherm constants
            eqilibration... currently does nothing
            algorithm... Choose between "Nelder-Mead" and "CMA-ES"
                the class will automatically create a minimisation problem based on the chosen algorithm
            
            tolf... the function value to cutoff, used only in Nelder-Mead
            tolx... the difference in x for the algorithm to cutoff on, used by both algorithms
        
        algortim = "CMA-ES raw" is available but should not be used as it performs worse 
        than the one used by the cma package
        """
          
        self.loss_function=electrokitty_loss(self.kin, self.species_information, self.cell_const
                                          ,self.isotherm, self.I_data,
                                          fit_Cdl=fit_Cdl, fit_Ru=fit_Ru, fit_gamamax=fit_gamamax,
                                          fit_A=fit_A, fit_iso=fit_iso)
        
        self.tells, self.gammaposition = self.loss_function.give_tells_gp()
        self.simulator.create_optimization_problem(self.tells, self.gammaposition)
        self.loss_function.update_ysim(self.simulator.calc_from_guess)
        
        if algorithm != "Nelder-Mead":
            lower_bound, upper_bound = self.loss_function.create_lower_upper_bounds(self.loss_function.guess, self.tells,
                                                                                      self.E_generated)
        else:
            lower_bound=None
            upper_bound=None
        
        self.optimizer=electrokitty_optimization_controller(self.loss_function.RMSE,
                                                            self.loss_function.guess,
                                                            algorithm=algorithm, 
                                                            tolf=tolf, tolx=tolx,
                                                            lb=lower_bound[:-1], ub=upper_bound[:-1])
        
        opt_params, self.fit_score=self.optimizer.fit_parameters()
        self.update_after_min(opt_params)
        
        current_time=str(datetime.now())
        current_time=current_time.replace(" ", "_")
        current_time=current_time.replace(":", "-")

        self.save(current_time)
        print()
        print("Finished Optimization and updated problem")
        
    def fit_harmonics(self,  base_freq, N_harmonics, w, fit_Cdl=False, fit_Ru=False, fit_gamamax=False,
                    fit_A=False, fit_iso=False, eqilibration=False, algorithm="Nelder-Mead",
                    tolf=10**-11, tolx=10**-11):
        """
        Function which tries to fit kinetic parameters and others to current harmonics generated from data
        
        Parameters:
            base_freq... frequency of the sine wave [Hz]
            N_harmonics... number of harmonics to use when fitting
            w... a list of vaules to use in the rectangular function to seperate the harmonics
            fit_Cdl... True to fit the double layer capacitance
            fit_Ru... True to fit uncompensated resistance
            fit_gamamax... True to fit the highest surface concentration in the species_information list
            fit_A... True to fit the electrode surface
            fit_iso... True to fit the whole list of isotherm constants
            eqilibration... currently does nothing
            algorithm... Choose between "Nelder-Mead" and "CMA-ES"
                the class will automatically create a minimisation problem based on the chosen algorithm
            
            tolf... the function value to cutoff, used only in Nelder-Mead
            tolx... the difference in x for the algorithm to cutoff on, used by both algorithms
        
        algortim = "CMA-ES raw" is available but should not be used as it performs worse 
        than the one used by the cma package
        """
      
        self.FFT_analyze_data(base_freq, N_harmonics, w)
        self.loss_function=electrokitty_loss(self.kin, self.species_information, self.cell_const
                                          ,self.isotherm, self.I_data,
                                          fit_Cdl=fit_Cdl, fit_Ru=fit_Ru, fit_gamamax=fit_gamamax,
                                          fit_A=fit_A, fit_iso=fit_iso)
        
        self.tells, self.gammaposition = self.loss_function.give_tells_gp()
        self.simulator.create_optimization_problem(self.tells, self.gammaposition)
        self.loss_function.update_ysim(self.simulator.calc_from_guess)
        
        self.loss_function.create_ACV_problem(base_freq, N_harmonics, self.I_har_data, self.t, w=w)
        
        if algorithm != "Nelder-Mead":
            lower_bound, upper_bound = self.loss_function.create_lower_upper_bounds(self.loss_function.guess, self.tells,
                                                                                      self.E_generated)
            lower_bound=lower_bound[:-1]
            upper_bound=upper_bound[:-1]
        else:
            lower_bound=None
            upper_bound=None
        
        self.optimizer=electrokitty_optimization_controller(self.loss_function.RMSE_har,
                                                            self.loss_function.guess,
                                                            algorithm=algorithm, 
                                                            tolf=tolf, tolx=tolx,
                                                            lb=lower_bound, ub=upper_bound)
        
        opt_params, self.fit_score=self.optimizer.fit_parameters()
        self.update_after_min(opt_params)
        self.FFT_analyze_sim(base_freq, N_harmonics, w)
        
        current_time=str(datetime.now())
        current_time=current_time.replace(" ", "_")
        current_time=current_time.replace(":", "-")

        self.save(current_time)
        print()
        print("Finished Optimization and updated problem")
    
    def update_after_min(self, optimal):
        """
        Function that updates class parameters after it finished fitting
        All parameters are overwritten
        """
        kine, cells, spinfo, isot = self.loss_function.unpack_fit_params(optimal, self.tells, self.gammaposition)
        
        self.kin=kine
        self.cell_const=cells
        self.species_information=spinfo
        self.isotherm=isot
        
        spectators = [np.ones(len(self.species_information[0])),np.ones(len(self.species_information[1]))]

        self.mechanism_list=self.Parser.Parse_mechanism()
        self.simulator.set_parameters(
                              cells, self.diffusion_const, isot, spectators, self.spatial_info, spinfo, kine, 
                              self.mechanism_list[0], self.mechanism_list[1], 
                              self.mechanism_list[2], self.mechanism_list[3], self.mechanism_list[4]
                              
                              )

        self.simulator.set_simulation_programm(self.t, self.E_generated)
        
        
        self.current = self.simulator.simulate()
        self.E_Corr = self.simulator.give_E_corr()
        self.surface_profile = self.simulator.give_surf_profile()
        self.concentration_profile = self.simulator.give_concentration_profile()

       
    def sample_parameter_distribution(self, n_samples=2000, burn_in_per=0.3, num_chains=1, multi_processing=False,
                                      fit_Cdl=False, fit_Ru=False, fit_gamamax=False,
                                                      fit_A=False, fit_iso=False, eqilibration=False, bounds=None):
        
        """
        Function which tries to fit kinetic parameters and others to data, using MCMC
        
        Currently multiprocessing does not fully work and should be avoided
        
        Parameters:
            n_samples... number of samples in a single chain
            burn_in_per... fraction of samples to be discared when updating the class
            num_chains... number of chains to be calculated
            multiprocessing... True in order to perform MCMC on multiple cores
            fit_Cdl... True to fit the double layer capacitance
            fit_Ru... True to fit uncompensated resistance
            fit_gamamax... True to fit the highest surface concentration in the species_information list
            fit_A... True to fit the electrode surface
            fit_iso... True to fit the whole list of isotherm constants
            eqilibration... currently does nothing
            bounds... a list containing a lists of the lower and upper bounds on the parameters,
                if None the class generates them itself, used as the prior distribution
            
        """
        
        self.multi_core_MCMC=multi_processing
        
        self.loss_function=electrokitty_loss(self.kin, self.species_information, self.cell_const
                                         ,self.isotherm, self.I_data,
                                         fit_Cdl=fit_Cdl, fit_Ru=fit_Ru, fit_gamamax=fit_gamamax,
                                         fit_A=fit_A, fit_iso=fit_iso)
        
        self.tells, self.gammaposition = self.loss_function.give_tells_gp()
        self.simulator.create_optimization_problem(self.tells, self.gammaposition)
        self.loss_function.update_ysim(self.simulator.calc_from_guess)
        if bounds==None:
            lower_bound, upper_bound = self.loss_function.create_lower_upper_bounds(self.loss_function.guess, self.tells,
                                                                                 self.E_generated)
            bounds=[lower_bound, upper_bound]
        self.MCMC_sampler=electrokitty_sampler(n_samples, burn_in_per, num_chains, 
                     multi_processing, bounds, self.I_data)
        
        self.MCMC_sampler.give_y_sim(self.simulator.calc_from_guess)
        
        chains=self.MCMC_sampler.start_sampler(self.loss_function.guess)
        
        
        self.chains=chains
        self.mean_chain=np.average(chains, axis=0)
        
        self.update_after_min(np.average(self.mean_chain[:int(burn_in_per*n_samples),:-1],axis=0))
        
        now=str(datetime.now)
        now.replace(" ", "_")
        now.replace(":", "-")
        self.save(now)
        
        print()
        print(f"Finished sampling after iterations: {n_samples}")
        print(f"Parameters estimated for burn in being: {100*burn_in_per}% of samples")
    
        
    ##################### Functions for treating and ploting data
    def FFT_analysis(self, f,N,w, current):
        """
        Function which generates the FT of the data and separates the harmonics
        
        Parameters:
            f... frequency of the AC_wave [Hz]
            N... number of harmonics to be extracted
            w... a list containing widths for a rectangular function when separating harmonics [Hz]
            current... the current array to be analysed, the class does this by itself, if using FFT_anlyse_data/sim
        """
        
        def rectangular(f,w0,w):
            return np.where(abs(f-w0)<=w,1,0)
        
        I_harmonics=[]
        dt=np.average(np.diff(self.t))
        freq=np.fft.fftfreq(self.t.shape[-1],d=dt)
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
        return sp, freq, I_harmonics
    
    def FFT_analyze_data(self, f, N, w):
        """
        Function that uses FT to analyse data and get harmonic currents, as well as the FT of the data
        
        Parameters:
            f... base frequency to use when seprateting harmonics
            N... number of harmonics to generate
            w... a list containing the widths for the rectangular function
        """
        self.sp, self.freq, self.I_har_data=self.FFT_analysis(f, N, w, self.I_data)
    
    def FFT_analyze_sim(self, f, N, w):
        """
        Function that uses FT to analyse simulation and get harmonic currents, as well as the FT of the data
        
        Parameters:
            f... base frequency to use when seprateting harmonics
            N... number of harmonics to generate
            w... a list containing the widths for the rectangular function
        """
        self.sp, self.freq, self.I_harmonics=self.FFT_analysis(f, N, w, self.current)
    
    def Harmonic_plots(self, plot_sim=False, plot_data=False,w=0,label="",Title=""):
        """
        Function to plot the generated harmonics. The output is N plots, each beeing a harmonic
        
        Parameters:
            plot_sim... True to plot the simulated harmonics
            plot_data... True to plot harmonics form data
            
            !one of these must be True
            
            w... choosing what kind of plot to make
                0 - i_har vs.t
                1 - i_har vs. E_dc
        """
        if plot_sim:
            I_hars=self.I_harmonics
        elif plot_data:
            I_hars=self.I_har_data    
        else:
            print("choose data or sim")
        
        def rectangular(f,w0,w):
            return np.where(abs(f-w0)<=w,1,0)
        
        if w==0:
            for i in range(len(self.I_harmonics)):
                plt.figure(str(i)+"_harmonic")
                plt.title(str(i)+"th harmonic "+Title)
                plt.plot(self.t, I_hars[i],label=label)
                plt.xlabel("time [s]")
                plt.ylabel("I [A]")
                # plt.savefig(str(i)+"_harmonic")
        
        elif w<0:
            print("Error in ploting function, w < 0")
        
        elif w==1:
            sp=np.fft.fft(self.E_generated)
            dt=self.t[1]-self.t[0]
            freq=np.fft.fftfreq(self.E_generated.shape[-1],d=dt)
            E_sim=np.fft.ifft(rectangular(freq,0,w)*sp).real
            for i in range(len(I_hars)):
                plt.figure(str(i)+"_harmonic")
                plt.title(str(i)+"th harmonic "+Title)
                plt.plot(E_sim,I_hars[i],label=label)
                plt.xlabel("$E_{DC}$ [V]")
                plt.ylabel("I [A]")
                # plt.savefig(str(i)+"_harmonic")
    
    def FT_plot(self,Title="",label=""):
        """
        Function plots the FT of either data of simulation, depending which was made last
        
        The plot is log10(abs(FT(i))) vs freq
        """
        plt.figure("Fourier_Transform")
        plt.title("Fourier Transform "+Title)
        plt.plot(self.freq,np.log10(np.abs(self.sp)),label=label)
        plt.xlabel("frequency [Hz]")
        plt.ylabel("$log_{10}$($I_{FT}$) [dB]")
        
    def Plot_simulation(self, Title="",label="",x_label=""):
        """
        Function plots the simulated current vs. potential used in the simulation
        """
        plt.figure("Base_signal")
        plt.title("Simulated signal "+Title)
        plt.plot(self.E_generated,self.current,label=label)
        plt.xlabel(x_label)
        plt.ylabel("I [A]")
        plt.show()
        
    def Plot_data(self, Title="",label="",x_label=""):
        """
        Function plots the data current vs. data potential
        """
        plt.figure("Base_signal")
        plt.title("Simulated signal "+Title)
        plt.plot(self.E_generated,self.I_data,label=label)
        plt.xlabel(x_label)
        plt.ylabel("I [A]")
        plt.show()
    
    def _update_conc_profile(self):
        """
        not to be used
        updates the consentration profile so that it can be plotted
        """
        waste_list=[]
        for i in range(self.number_of_diss_spec):
            waste_list.append(self.concentration_profile[:,i:-self.number_of_diss_spec+1:self.number_of_diss_spec])
        return waste_list
    
    def Plot_concentration_profile(self, species_num=0):
        """
        Function plots the concentration profile as a 3D plot
        
        The plot is as (t, x, c)
        
        Parameters:
            species_num... the list index of the species to be plotted
            
            to see which index a species has, you can call <class_name>.mechanism_list[0][1]
        """
        
        def _find_gama(dx,xmax,nx):
            # bisection method for finding gama
            # used in determining the exponential spatial grid
            a=1
            b=2
            for it in range(0,50):
                gama=(a+b)/2
                f=dx*(gama**nx-1)/(gama-1)-xmax
                if f<=0:
                    a=gama
                else:
                    b=gama
                if abs(b-a)<=10**-8:
                    break
            gama=(a+b)/2
            if gama>2:
                print("bad gama value")
                sys.exit()
            return gama
        
        def calc_x(tmax,D,fraction,nx):
            # Given the simulation time, f, the maximum diffusion coefficient, the initial dx
            # and the lenghth of spatial direction
            # evaluates a one dimensional grid to be used in simulation
            # fraction is given as dx/xmax
            xmax=6*np.sqrt(tmax*D)
            dx=fraction*xmax
            gama=_find_gama(dx, xmax, nx)
            N=np.arange(nx+2)
            x=dx*(gama**N-1)/(gama-1)
            return x
        
        
        try:
            self.concentration_profile=self._update_conc_profile()
        except:
            pass
        self.x=calc_x(self.t[-1], max(self.diffusion_const), self.spatial_info[0], self.spatial_info[1])
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

        X, Y = np.meshgrid(self.x[:-2]*10**3, self.t)
        Z=self.concentration_profile[species_num][:,:len(self.x[:-2])]
        
        surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)
        
        ax.set_zlim(0,np.max(self.concentration_profile[species_num]))
        
        ax.set_xlabel("Electrode_Distance [mm]")
        ax.set_ylabel("Time [s]")
        ax.set_zlabel("concentration [mol/m^3]")
        ax.set_title("concentration profile for: "+self.mechanism_list[0][1][species_num])
        
        plt.show()
        
    def Plot_Adsorbed_species(self,Title="",label="",x_label="",mode=0):
        """
        Function plots the surface concentration profile
        The plot is as: gama_i vs. potential in simulation
        
        The function plots all surface concentrations
        """
        plt.figure("Adsorbed_species")
        plt.title("Adsorbed_species "+Title)
        if mode==0:
            plt.plot(self.E_generated,self.surface_profile,label=label)
            plt.xlabel(x_label)
            plt.ylabel("theta")
            plt.legend(self.mechanism_list[0][0])
            plt.show()
        else:
            plt.plot(self.t,self.surface_profile,label=label)
            plt.xlabel(x_label)
            plt.ylabel("theta")
            plt.show()
            