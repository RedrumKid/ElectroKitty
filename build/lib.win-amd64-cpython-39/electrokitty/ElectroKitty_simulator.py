# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 10:47:28 2024

@author: ozbejv
"""

import sys
import itertools
import numpy as np
import scipy.optimize as sciop
from. ElectroKitty_parameter_distributions import gaussian_distribution
from cpp_ekitty_simulator import cpp_ekitty_simulator

class electrokitty_simulator:
    
    def __init__(self):
        self.i_mean = None
        self.sim_cell_const = None
        self.sim_diffusion_const = None
        self.sim_isotherm = None
        self.sim_spectators = None
        self.sim_spatial_info = None
        self.sim_species_information = None
        self.sim_kin = None
        self.simulate_with_dispersion = False
        self.dispersed_species_information = None
        self.dispersed_cell_const = None
        self.dispersed_isotherm = None
        self.dispersed_kin = None

        self.tells = None
        self.gamapos = None
        
        self.sim_mechanism_list = None
        self.save_kin = None
        self.save_spec_info = None
        self.save_cell_cons = None
        self.save_isotherm = None
        
        self._simulator = None
    
    def check_type(self, param):
        if type(param) is float or type(param) is int:
            return True
        else:
            return False
    
    def check_for_params(self, some_list):
        scratch_list = []
        dispersion_check = False
        for ind in range(len(some_list)):
            if self.check_type(some_list[ind]):
                scratch_list.append(True)
            else: 
                scratch_list.append(False)
                dispersion_check = True
        return scratch_list, dispersion_check
    
    def create_disp_lists(self):
        self.dispersed_cell_const, check = self.check_for_params(self.sim_cell_const)
        if check:
                self.simulate_with_dispersion = check

        self.dispersed_isotherm, check = self.check_for_params(self.sim_isotherm)
        if check:
                self.simulate_with_dispersion = check

        self.dispersed_species_information = []
        for spec_in in self.sim_species_information:
            in_between, check = self.check_for_params(spec_in)
            self.dispersed_species_information.append(in_between)
            if check:
                self.simulate_with_dispersion = check
        
        self.dispersed_kin = []
        for ki in self.sim_kin:
            in_between, check = self.check_for_params(ki)
            self.dispersed_kin.append(in_between)
            if check:
                self.simulate_with_dispersion = check
            
    
    def give_simulation_constants(self, kins, cell_consts, 
                          Diffusion_consts, isotherms ,Spatial_infos , 
                          Species_informations, spectatorss=False):
        
        self.sim_cell_const=cell_consts
        self.sim_diffusion_const=Diffusion_consts
        self.sim_isotherm=isotherms
        self.sim_spectators=spectatorss
        self.sim_spatial_info=Spatial_infos
        self.sim_species_information=Species_informations
        self.sim_kin=kins
        spectators = [np.ones(len(Species_informations[0])),np.ones(len(Species_informations[1]))]
        self.sim_spectators = spectators

        self.create_disp_lists()
    
    def give_mechanism_list(self, mechanism_list):
        self.sim_mechanism_list = mechanism_list
    
    def give_simulation_program(self, ts, E_gens):
        self.t = ts
        self.E_gen = E_gens

    def give_saved_cons(self):
        return self.save_kin, self.save_spec_info, self.save_cell_cons, self.save_isotherm

    def create_weights(self, fun, x):
        ws = []
        points = []
        for i in range(1, len(x)):
            ws.append((x[i]-x[i-1])*fun((x[i]+x[i-1])/2))
            points.append((x[i]+x[i-1])/2)
        return np.array(ws), np.array(points)

    def create_integration_scale(self, xmin, xmax, N):
        return np.linspace(xmin, xmax, N+1)

    def create_dist_simulation_list(self):
        # 0-kin; 1-spec_info; 2-cell_const; 3-iso
        simulation_list = []
        ws = []
        xs = []
        # check kin
        for i in range(len(self.dispersed_kin)):
            for j in range(len(self.dispersed_kin[i])):
                if self.dispersed_kin[i][j] == False:
                    w, x = self.create_weights(self.sim_kin[i][j][0], 
                                               self.create_integration_scale(self.sim_kin[i][j][2], self.sim_kin[i][j][3], self.sim_kin[i][j][1]))
                    if self.sim_kin[i][j][4] == "log":
                        x = np.exp(x)
                    simulation_list.append([0, [i, j]])
                    ws.append(w)
                    xs.append(x)
        
        # check species_information
        for i in range(len(self.dispersed_species_information)):
            for j in range(len(self.dispersed_species_information[i])):
                if self.dispersed_species_information[i][j] == False:
                    w, x = self.create_weights(self.sim_species_information[i][j][0], 
                                               self.create_integration_scale(self.sim_species_information[i][j][2], 
                                                                             self.sim_species_information[i][j][3], self.sim_species_information[i][j][1]))
                    if self.sim_species_information[i][j][4] == "log":
                        x = np.exp(x)
                    simulation_list.append([1, [i, j]])
                    ws.append(w)
                    xs.append(x)

        # check cell_const
        for i in range(len(self.dispersed_cell_const)):
            if self.dispersed_cell_const[i] == False:
                w, x = self.create_weights(self.sim_cell_const[i][0], self.create_integration_scale(self.sim_cell_const[i][2],
                                                                                                self.sim_cell_const[i][3], self.sim_cell_const[i][1]))
                if self.sim_cell_const[i][4] == "log":
                    x = np.exp(x)
                simulation_list.append([2, [i, 0]])
                ws.append(w)
                xs.append(x)
        
        # check isotherm
        for i in range(len(self.dispersed_isotherm)):
            if self.dispersed_isotherm[i] == False:
                w, x = self.create_weights(self.sim_isotherm[i][0], self.create_integration_scale(self.sim_isotherm[i][2],
                                                                                              self.sim_isotherm[i][3], self.sim_isotherm[i][1]))
                if self.sim_isotherm[i][4] == "log":
                    x = np.exp(x)
                simulation_list.append([3, [i, 0]])
                ws.append(w)
                xs.append(x)
        return simulation_list, ws, xs

    def simulate_dispersion(self):
        simulation_list, ws, xs = self.create_dist_simulation_list()

        mean_i = 0
        mean_E_corr = 0
        mean_adsorbed_spec = 0
        mean_conc_prof = 0

        kin = list(self.sim_kin)
        cell_const = list(self.sim_cell_const)
        species_info = list(self.sim_species_information)
        iso = list(self.sim_isotherm)
        
        for combinations in zip(itertools.product(*ws), itertools.product(*xs)):
            product = 1
            for w in combinations[0]:
                product *= w

            for ind in range(len(combinations[1])):
                if simulation_list[ind][0] == 0:
                    kin[simulation_list[ind][1][0]][simulation_list[ind][1][1]] = combinations[1][ind]
                elif simulation_list[ind][0] == 1:
                    species_info[simulation_list[ind][1][0]][simulation_list[ind][1][1]] = combinations[1][ind]
                elif simulation_list[ind][0] == 2:
                    cell_const[simulation_list[ind][1][0]] = combinations[1][ind]
                elif simulation_list[ind][0] == 3:
                    iso[simulation_list[ind][1][0]] = combinations[1][ind]
                
            self.simulator = cpp_ekitty_simulator()
            self.simulator.set_parameters(
                                cell_const, self.sim_diffusion_const, iso, self.sim_spectators, 
                                self.sim_spatial_info, species_info, kin, 
                                self.sim_mechanism_list[0], self.sim_mechanism_list[1], 
                                self.sim_mechanism_list[2], self.sim_mechanism_list[3], self.sim_mechanism_list[4]
                                )

            self.simulator.set_simulation_programm(self.t, self.E_gen)
            
            current = self.simulator.simulate()
            E_Corr = self.simulator.give_E_corr()
            surface_profile = self.simulator.give_surf_profile()
            concentration_profile = self.simulator.give_concentration_profile()

            mean_i += product*current
            mean_E_corr += product*E_Corr
            mean_adsorbed_spec += product*surface_profile
            mean_conc_prof += product*concentration_profile

        return mean_i, mean_E_corr, mean_adsorbed_spec, mean_conc_prof
    
    def import_for_fitting(self, tells, gama_position):
        self.tells = tells
        self.gamapos = gama_position

    def calc_from_guess(self, guess):
        self.sim_kin, self.sim_cell_const, self.sim_species_information, self.sim_isotherm = self.unpack_fit_params(guess, self.tells, self.gamapos,
                                                                                                    self.sim_kin, self.sim_species_information, 
                                                                                                    self.sim_cell_const, self.sim_isotherm)
        
        i_sim, e, a, ds = self.simulate()
        return i_sim
    
    def unpack_fit_params(self, guess, tells, gamma_position, kins, species_informations, cell_consts, isotherms):
        """
        Function takes the guess, tells and gammma_position to reconstruct the lists for the simulator
        """
        guess=guess.tolist()
        kinetics=list(kins)
        cell_params=list(cell_consts)
        spec_info=list(species_informations)
        iso = list(isotherms)
        
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


    def simulate(self):

        if self.simulate_with_dispersion:
            current, E_Corr, surface_profile, concentration_profile = self.simulate_dispersion()

        else: 
            self.simulator = cpp_ekitty_simulator()
            self.simulator.set_parameters(
                                self.sim_cell_const, self.sim_diffusion_const, self.sim_isotherm, self.sim_spectators, 
                                self.sim_spatial_info, self.sim_species_information, self.sim_kin, 
                                self.sim_mechanism_list[0], self.sim_mechanism_list[1], 
                                self.sim_mechanism_list[2], self.sim_mechanism_list[3], self.sim_mechanism_list[4]
                                )

            self.simulator.set_simulation_programm(self.t, self.E_gen)
            
            current = self.simulator.simulate()
            E_Corr = self.simulator.give_E_corr()
            surface_profile = self.simulator.give_surf_profile()
            concentration_profile = self.simulator.give_concentration_profile()
        
        return current, E_Corr, surface_profile, concentration_profile

class python_electrokitty_simulator:
    """
    Python version of the simulator
    ElectroKitty uses a c++ implementation of this code
    mostly here for refrencing and testing
    """
    def __init__(self):
        self.F=96485
        self.R=8.314
        self.t=None
        self.E_generated=None
        self.current=None
        self.concentration_profile=None
        self.surface_profile=None
        
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
        self.tells=None
        self.gammaposition=None
    
    def update_parameters(self, mechanism_list, kin, cell_const, 
                          Diffusion_const, isotherm,Spatial_info, 
                          Species_information, spectators=False):
        
        self.cell_const=cell_const
        self.diffusion_const=Diffusion_const
        self.isotherm=isotherm
        self.spectators=spectators
        self.spatial_info=Spatial_info
        self.species_information=Species_information
        self.kin=kin
        
        self.mechanism_list=mechanism_list
    def give_sim_program(self, E, t):
        self.E_generated=E
        self.t=t
    
    ############################## Functions for precalc and simulator
    
    def _get_kinetic_constants(self, k_vector, kinetic_types):
        # A function wich checks the number of constants for reversible or irreversible steps and makes the lists correct (assigns a zero for the back step in irreversible steps)
        for i in range(len(k_vector)):
            if kinetic_types[i]==0 and len(k_vector[i])!=2:
                print("Error in number of constants: Reversible step, assigned 1 constants, requres 2")
                sys.exit()
                
            elif kinetic_types[i]==1 and len(k_vector[i])==1:
                test=[0]
                test.append(k_vector[i][0])
                k_vector[i]=test
                
            elif kinetic_types[i]==2 and len(k_vector[i])==1:
                k_vector[i].append(0)
                
            elif kinetic_types[i]==1 or kinetic_types[i]==2 and len(k_vector[i])!=1:
                print("Error in number of constants: Irreversible step, assigned more constants, requres 1")
                sys.exit()
        
        return k_vector
    
    def iterate_Over_conc(self, step, c, term, isotherm):
        
        for i in step:
            term*=c[i]*np.exp(-isotherm[i]*c[i])
        return term

    def update_K_Matrix(self, k_matrix, term_f, term_b, step):
        check=[]
        for i in step:
            if i not in check:
                check.append(i)
                k_matrix[i]+=-term_f
                k_matrix[i]+=term_b
        return k_matrix
     
    def calc_kinetics(self, reac_type,c,index,kinetic_const, isotherm):
        # A function that given the reaction type 0-ads, 1- bulk 
        # the relevant concentrations, indexes connecting the c and constants, 
        # evaluates the reaction forward and backward kinetic rate
        k_matrix=np.zeros(len(c))
    
        for i in range(len(index[reac_type])):
            
            step=index[reac_type][i]
            constants=kinetic_const[i]
            
            forward_step=constants[0]
            backward_step=constants[1]
            
            forward_step=self.iterate_Over_conc(step[0], c, forward_step, isotherm)
            backward_step=self.iterate_Over_conc(step[1], c, backward_step, isotherm)
            
            k_matrix=self.update_K_Matrix(k_matrix, forward_step, backward_step, step[0])
            k_matrix=self.update_K_Matrix(k_matrix, -forward_step, -backward_step, step[1])
                    
        return k_matrix

    def _calc_EC_kinetics(self, reac_type, c, index, kinetic_const, E, isotherm):
        # A function that given the reaction type 0-ads, 1- bulk 
        # the relevant concentrations, indexes connecting the c and constants, 
        # evaluates the reaction forward and backward electrochemical kinetic rate at the boundary
        k_matrix=np.zeros(len(c))
    
        for i in range(len(index[reac_type])):
            
            step=index[reac_type][i]
            constants=kinetic_const[i]
            
            forward_step=constants[0](E)
            backward_step=constants[1](E)
            
            forward_step=self.iterate_Over_conc(step[0], c, forward_step, isotherm)
            backward_step=self.iterate_Over_conc(step[1], c, backward_step, isotherm)
            
            k_matrix=self.update_K_Matrix(k_matrix, forward_step, backward_step, step[0])
            k_matrix=self.update_K_Matrix(k_matrix, -forward_step, -backward_step, step[1])
                    
        return k_matrix

    def _calc_current(self, reac_type, c, index, kinetic_const, E, isotherm):
        # A function that given the reaction type 0-ads, 1- bulk 
        # the relevant concentrations, indexes connecting the c and constants, 
        # evaluates the reaction forward and backward electrochemical current
        # the output must be multiplied with n*F*A
        current=0
    
        for i in range(len(index[reac_type])):
            step=index[reac_type][i]
            constants=kinetic_const[i]
            forward_step=constants[0](E)
            backward_step=constants[1](E)
            forward_step=self.iterate_Over_conc(step[0], c, forward_step, isotherm)
            backward_step=self.iterate_Over_conc(step[1], c, backward_step, isotherm)
            current+=-forward_step+backward_step
        return current
    
    def _find_gama(self, dx,xmax,nx):
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
    
    def _Fornberg_weights(self, z,x,n,m):
    # From Bengt Fornbergs (1998) SIAM Review paper.
    #  	Input Parameters
    #	z location where approximations are to be accurate,
    #	x(0:nd) grid point locations, found in x(0:n)
    #	n one less than total number of grid points; n must
    #	not exceed the parameter nd below,
    #	nd dimension of x- and c-arrays in calling program
    #	x(0:nd) and c(0:nd,0:m), respectively,
    #	m highest derivative for which weights are sought,
    #	Output Parameter
    #	c(0:nd,0:m) weights at grid locations x(0:n) for derivatives
    #	of order 0:m, found in c(0:n,0:m)
    #      	dimension x(0:nd),c(0:nd,0:m)
        
        c=np.zeros((n+1,m+1))
        c1=1
        c4=x[0]-z
        
        c[0,0]=1
        
        for i in range(1,n):
            mn=min([i,m])
            c2=1
            c5=c4
            c4=x[i]-z
            for j in range(0,i):
                c3=x[i]-x[j]
                c2=c3*c2
                
                if j==i-1:
                    for k in range(mn,0,-1):
                        c[i,k]=c1*(k*c[i-1,k-1]-c5*c[i-1,k])/c2
                    c[i,0]=-c1*c5*c[i-1,0]/c2
                
                for k in range(mn,0,-1):
                    c[j,k]=(c4*c[j,k]-k*c[j,k-1])/c3
                
                c[j,0]=c4*c[j,0]/c3
            
            c1=c2
        
        return c
    
    def _Space_ranges(self, tmax,f,D,fraction,nx):
        # Given the simulation time, f, the maximum diffusion coefficient, the initial dx
        # and the lenghth of spatial direction
        # evaluates a one dimensional grid to be used in simulation
        # fraction is given as dx/xmax
        xmax=6*np.sqrt(tmax*D)
        dx=fraction*xmax
        gama=self._find_gama(dx, xmax, nx)
        N=np.arange(nx+2)
        self.x=dx*(gama**N-1)/(gama-1)
        return self.x
    
    def _calc_main_coef(self, x,dt,D,nx,B):
        # calculate alfas and a's used in simulation
        # calculated with given spatial direction x
        # the weights are given via the method of finite difference implicit method
        # B is to be implemented
        a1=[]
        a2=[]
        a3=[]
        a4=[]
        
        for i in range(1,nx):
            
            weights=self._Fornberg_weights(x[i],x[i-1:i+3],4,2)
            
            alfa1d=weights[0,2]
            alfa2d=weights[1,2]
            alfa3d=weights[2,2]
            alfa4d=weights[3,2]
            
            alfa1v=-(B*x[i]**2)*weights[0,1]
            alfa2v=-(B*x[i]**2)*weights[1,1]
            alfa3v=-(B*x[i]**2)*weights[2,1]
            alfa4v=-(B*x[i]**2)*weights[3,1]
            
            a1.append((-alfa1d*D-alfa1v)*dt)
            a2.append((-alfa2d*D-alfa2v)*dt+1)
            a3.append((-alfa3d*D-alfa3v)*dt)
            a4.append((-alfa4d*D-alfa4v)*dt)
        
        return np.array([np.array(a1),np.array(a2),np.array(a3),np.array(a4)])

    def _calc_boundary_condition(self, x,i,D,nx,B):
        # A function for evaluation of the flux boundary condition, at either boundary
        # i should be 0 or -1, 0 for the electrode, -1 for the bulk limit
        # B is used in case of rotation (to be implemented)
        a1=[]
        a2=[]
        a3=[]
        
        if i==0:
            weights=self._Fornberg_weights(x[i],x[i:i+3],3,1)
        elif i==-1:
            weights=self._Fornberg_weights(x[i],x[i-2:],3,1)
        else:
            print("Boundary Error: boundary flux indexed incorrectly")
            
        alfa1=weights[0,1]-(B*x[i]**2)
        alfa2=weights[1,1]-(B*x[i]**2)
        alfa3=weights[2,1]-(B*x[i]**2)
        
        a1.append(-alfa1*D)
        a2.append(-alfa2*D)
        a3.append(-alfa3*D)
        
        return np.array([np.array(a1),np.array(a2),np.array(a3)])
    
    def _Butler_volmer_kinetics(self, alpha, k0, E0, f, el_num):
        # A function for evaluating the butler-volmer kinetics 
        # it transforms the given constants into function to be evaluated during simulation
        return [lambda E: el_num*k0*np.exp(-alpha*el_num*f*(E-E0)), lambda E:el_num*k0*np.exp((1-alpha)*el_num*f*(E-E0))]
    
    def _get_EC_kinetic_constants(self, k_vector, kinetic_types, f, num_el):
        # A function for getting BV kinetics at the boundary condition
        # in case of irreversible kinetics the function is a zero function
        for i in range(len(k_vector)):
            k_vector[i]=self._Butler_volmer_kinetics(k_vector[i][0], k_vector[i][1], k_vector[i][2], f, num_el[i])
            if kinetic_types[i]==1:
                k_vector[i][0]=lambda E: 0
            elif kinetic_types[i]==2: 
                k_vector[i][1]=lambda E: 0
        return k_vector
    
    def _time_step(self, c, a, cp, nx, dt, n1, n, bound1, bound2, pnom, constants, index, F, delta, isotherm_constants, null, spectator):
        # A function for evaluating the time step
        # given the guess, the weights, previous iteration, number of x points,
        # dt, number of ads spec, number of bulk spec, boundary at the electrode
        # boundary at the bulk limit, the program value of potential, a list of constants ordered:
            # ads, bulk, ec, cell
        # the index of how are kinetics manipulated, faraday constant, and the derivative of the potential
        # evaluates the nonlinear set of equations to be solved at each time step
        Ru,Cdl,A=constants[-1][1:]
        p=c[-2]
    
        gc=c[-1]
        gcp=cp[-1]
        
        theta=c[:n1]
        thetap=cp[:n1]
        
        c=c[n1:-2]
        cp=cp[n1:-2]
    
        c=c.reshape((nx+2,n))
        cp=cp.reshape((nx+2,n))
    
        f=np.zeros(n1+(n)*(nx+2))
    
        bound_kinetics=(self.calc_kinetics(0, np.append(theta, c[0,:]), index, constants[0], isotherm_constants)
                        + self._calc_EC_kinetics(2,np.append(theta, c[0,:]), index, constants[2], p, isotherm_constants))
    
        f[:n1]=theta-thetap-dt*bound_kinetics[:n1]*spectator[:n1]
        
        f[n1:n1+n]=np.sum(bound1[:,0,:]*c[0:3,:],axis=0)-bound_kinetics[n1:]*spectator[n1:n1+n]
    
        if n!=0:
            for xx in range(1,nx):
                f[n1+n*xx:n1+n*xx+n]=(np.sum(a[:,xx-1,:]*c[xx-1:xx+3,:],axis=0)
                                      -dt*self.calc_kinetics(1, c[xx,:], index, constants[1], null)-cp[xx,:])
            f[-2*n:-n]=(c[-2,:])-bound2
            f[-n:]=(c[-2,:]-c[-1,:])
        else:
            pass
            
        ga=F*A*self._calc_current(2, np.append(theta, c[0,:]), index, constants[2], p, isotherm_constants)
        
        f9=(1+Ru*Cdl/dt)*gc-Cdl*delta-Ru*Cdl*(gcp)/dt
        f10=pnom-p-Ru*ga-Ru*gc
        f=np.append(f,np.array([f9,f10]))
        return f
    
    def _eqilibration_step(self, c, a, cp, nx, dt, n1, n, bound1, bound2, pnom, constants, index, F, delta, isotherm_constants, null, spectator):
        # A function for evaluating the time step
        # given the guess, the weights, previous iteration, number of x points,
        # dt, number of ads spec, number of bulk spec, boundary at the electrode
        # boundary at the bulk limit, the program value of potential, a list of constants ordered:
            # ads, bulk, ec, cell
        # the index of how are kinetics manipulated, faraday constant, and the derivative of the potential
        # evaluates the nonlinear set of equations to be solved at each time step
        
        Ru,Cdl,A=constants[-1][1:]
        p=c[-2]
    
        gc=c[-1]
        gcp=cp[-1]
        
        theta=c[:n1]
        
        c=c[n1:-2]
        cp=cp[n1:-2]
    
        c=c.reshape((nx+2,n))
        cp=cp.reshape((nx+2,n))
    
        f=np.zeros(n1+(n)*(nx+2))
    
        bound_kinetics=(self.calc_kinetics(0, np.append(theta, c[0,:]), index, constants[0], isotherm_constants)
                        + self._calc_EC_kinetics(2,np.append(theta, c[0,:]), index, constants[2], p, isotherm_constants))
    
        f[:n1]=bound_kinetics[:n1]*spectator[:n1]
        
        f[n1:n1+n]=bound_kinetics[n1:]*spectator[n1:n1+n]
    
        if n!=0:
            for xx in range(1,nx):
                f[n1+n*xx:n1+n*xx+n]=(np.sum(a[:,xx-1,:]*c[xx-1:xx+3,:],axis=0)/dt
                                      -self.calc_kinetics(1, c[xx,:], index, constants[1], null)-cp[xx,:])
            f[-2*n:-n]=(c[-2,:])-bound2
            f[-n:]=(c[-2,:]-c[-1,:])
        else:
            pass
            
        ga=F*A*self._calc_current(2, np.append(theta, c[0,:]), index, constants[2], p, isotherm_constants)
        
        f9=(1+Ru*Cdl/dt)*gc-Cdl*delta-Ru*Cdl*(gcp)/dt
        f10=pnom-p-Ru*ga-Ru*gc
        f=np.append(f,np.array([f9,f10]))
        return f
    
    def _create_const_list(self,indexs, const):
        c=[]
        for i in indexs:
           c.append(const[i])
        return c
    
    def simulator_Main_loop(self, mechanism_list, kin_const ,Constants, Spatial_info, Time, Species_information, Potential_program, eqilibration=True):
        # The main simulation function
        # Given the mechanism string given as 
            # C: or E: sum: f1 = or - sum b1 \n ...
        # A list of constants: list of lists for ads, bulk, ec, cell and diffusion
            # cell constants are supposed to be given as temperature, resistance, capacitance, electrode area 
        # Spatial info is a list with the fraction and number of x points, rest is evaluated by default
        # Time is requiered to be evenly spaced and is given as a numpy array
        # Species information is a list of two lists:
            # first contains the initial condition for adsorbed species given in gamas (moles per surface)
            # second is a list of functions to evaluate the initial condition of the concentration profile at t=0
        # Potential program is as a numpy array for wich the current is then updated
        
        # The function returns 3 arrays in given order: potential, current, time
        
        self.t=Time
        spec, index, types, r_ind, num_el=mechanism_list
        
        n=len(spec[1])
        n1=len(spec[0])
        self.number_of_surf_conf=n1
        self.number_of_diss_spec=n
        cell_const, Diffusion_const, isotherm_constants, spectator = Constants
        
        if spectator == False:
            spectator=np.ones(n+n1)
        else:
            spectator=np.array(spectator[0]+spectator[1])
        
        Diffusion_const=np.array(Diffusion_const)
        
        # isotherm_constants=isotherm_constants+n*[0]
        null=np.zeros(n)
        if n1>0:
            isotherm_constants=np.array(isotherm_constants)/max(Species_information[0])
        else:
            isotherm_constants=np.array(isotherm_constants)
        isotherm_constants=np.append(isotherm_constants,np.zeros(n))
        
        ads_const=self._create_const_list(r_ind[0], kin_const)
        bulk_const=self._create_const_list(r_ind[1], kin_const)
        EC_const=self._create_const_list(r_ind[2], kin_const)
        
        T,Ru,Cdl,A=cell_const
        f=self.F/self.R/T

        ads_const=self._get_kinetic_constants(ads_const, types[0])
        bulk_const=self._get_kinetic_constants(bulk_const, types[1])
        EC_const=self._get_EC_kinetic_constants(EC_const, types[2], f, num_el)

        dt=np.average(np.diff(Time))
        
        if len(Species_information[1])>0:
            self.x=self._Space_ranges(Time[-1], f, max(Diffusion_const), Spatial_info[0], Spatial_info[1])
        else:
            self.x=self._Space_ranges(Time[-1], f, 1, Spatial_info[0], Spatial_info[1])
        
        viscosity, ni = Spatial_info[2:]
        velocity_c=-0.51/np.sqrt(viscosity)*(2*np.pi*ni)**1.5
        
        a=self._calc_main_coef(self.x, dt, Diffusion_const, len(self.x)-2, velocity_c)
        
        theta=np.array(Species_information[0])
        
        c=np.zeros((len(self.x),n))
        for i in range(len(spec[1])):
            c[:,i]=Species_information[1][i]*np.ones(len(self.x))
        bound2=c[-1,:] 
        bound1=self._calc_boundary_condition(self.x, 0, Diffusion_const, 3, velocity_c)
        
        
        c=c.reshape((1,n*len(self.x)))[0,:]
        c=np.append(theta,c)
        
        delta_E=np.diff(Potential_program)/dt
        c=np.append(c,np.array([Potential_program[0],0]))
        
        constants=[ads_const, bulk_const, EC_const, cell_const]
        
        if eqilibration==True:
            # Preqilibration
            cp=c
            cp[-2]=Potential_program[0]
            
            res=sciop.root(self._eqilibration_step, cp, args=(
                a,cp,len(self.x)-2, dt , n1, n, 
                bound1, bound2, Potential_program[0], 
                constants, index, self.F, delta_E[0], isotherm_constants, null, spectator),tol=10**-28)
        
            c=res.x
            
            current=[]
            cap_cur=[]
            ps=[]
        
        else:
            current=[]
            cap_cur=[]
            ps=[]
        surface_profile=[]
        concentration_profile=[]
        for tt in range(0,len(Time)):
            cp=c
            cp[-2]=Potential_program[tt]
            res=sciop.root(self._time_step, cp, args=(
                a,cp,len(self.x)-2, dt , n1, n, 
                bound1, bound2, Potential_program[tt], 
                constants, index, self.F, delta_E[tt-1], isotherm_constants, null, spectator),tol=10**-28,method="hybr")
    
            c=res.x
            current.append(self.F*A*self._calc_current(2, c[:n1+n], index, EC_const, c[-2], isotherm_constants))
            cap_cur.append(c[-1])
            ps.append(c[-2])
            surface_profile.append(c[:n1])
            concentration_profile.append(c[n1:-2])
            
        ps=np.array(ps)
        current=np.array(current)
        cap_cur=np.array(cap_cur)
        current=current+cap_cur
        surface_profile=np.array(surface_profile)
        concentration_profile=np.array(concentration_profile)
        return ps, current, surface_profile, concentration_profile
    
    
    def simulate(self, eqilib=False):
        
        self.E_Corr, self.current, self.surface_profile, self.concentration_profile = self.simulator_Main_loop(
            self.mechanism_list, 
            self.kin, 
            [self.cell_const,
            self.diffusion_const,
            self.isotherm,
            self.spectators], 
            self.spatial_info, 
            self.t, 
            self.species_information, 
            self.E_generated, eqilibration=eqilib)
        
        return self.E_Corr, self.current, self.surface_profile, self.concentration_profile
