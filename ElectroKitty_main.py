# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 14:43:51 2024

@author: ozbejvodeb
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.optimize as sciop
import scipy.signal as scisi

################# Functions for parsing mechanism

# A script of functions for parsing a mechanism string
# Zero library dependancies requierd, for Python 3.x 
# Mechanistic steps must be separated by \n, the species in forward and backward reaction are
# separated by +, with = stating that the reaction is in equilibrium, and - saying 
# the reaction is one-directional
# Besides the characters, \n, +, =, -, any character or sequence of characters can be used,
# With the caveat that all spaces are discraded

def update_species_lists(line, bulk_spec, ads_spec, string):
    # A function for correctly identifying specis and sorting them into their respective lists
    for reac_fb in line.split(string):
        for s in reac_fb.split("+"):
            if "*" == s[-1]:
                if s not in ads_spec:
                    ads_spec.append(s)
            else:
                if s not in bulk_spec:
                    bulk_spec.append(s)
    return bulk_spec, ads_spec

def find_string_index(string, lis):
    # Find index of an element in a list
    i=0
    for element in lis:
        if element == string:
            return i
        else:
            i+=1

def get_index(string, master_list):
    # A function for finding and correctly indexing the species
    if string in master_list[0]:
        return find_string_index(string, master_list[0])
    elif string in master_list[1]:
        return find_string_index(string, master_list[1])+len(master_list[0])
    else:
        print("Mechanism error: Cannot find species in species list")
        return False

def index_list(line, species_list):
    # Finds the species indexs in the species list
    f_index=[]
    
    for spec in line.split("+"):
        f_index.append(get_index(spec, species_list))
        
    return f_index

def update_index_lists(ads_index, bulk_index , ec_index, f,b, line, reaction_types_ads, reaction_types_ec, reaction_types_bulk, Type, len_ads):
    
    # A function, which takes the index lists and updates them accordingly,
    # separating the reactions into either into adsorbtion, bulk or electrochemical
    
    if "*" in line and line[:2]=="C:":
        ads_index.append([f,b])
        reaction_types_ads.append(Type)
    elif line[:2]=="E:":
        ec_index.append([f,b])
        reaction_types_ec.append(Type)
    else:
        for i in range(len(f)):
            f[i]-=len_ads
        for i in range(len(b)):
            b[i]-=len_ads
        bulk_index.append([f,b])
        reaction_types_bulk.append(Type)
    return ads_index, bulk_index, ec_index, reaction_types_ads, reaction_types_ec, reaction_types_bulk

def Parse_mechanism(string):
    
    # The main function for analysing the mechanism string.
    # # The function takes the string, with each line (separated by \n),
    # each line must be declared as either a C or E mechanism, followed by :,
    # * at the end of the species, assignes the species as adsorbed,
    # and gives out a list of species, first list beiing the adsorbed, second the soulution,
    # and three lists of indecies to connect the species via differential equetions.
    # First are the adsobed species list, second is the bulk and finally the ec connections
    
    
    
    string=string.replace(" ","")
    a=string.split("\n")

    bulk_spec=[]
    ads_spec=[]

    bulk_index=[]
    ads_index=[]
    ec_index=[]
    
    reaction_types_ads=[]
    reaction_types_ec=[]
    reaction_types_bulk=[]
    
    for line in a:
        if "=" in line:
            bulk_spec, ads_spec = update_species_lists(line[2:], bulk_spec, ads_spec, "=")
        elif "-" in line:
            bulk_spec, ads_spec = update_species_lists(line[2:], bulk_spec, ads_spec, "-")
        else:
            print("Mechanism Error: Wrong mechanism separator")

    species_list=[ads_spec,bulk_spec]
    
    for line in a:
        if "=" in line:
            f,b=line[2:].split("=")
            f_ind=index_list(f, species_list)
            b_ind=index_list(b, species_list)
            ads_index, bulk_index, ec_index, reaction_types_ads,reaction_types_ec,reaction_types_bulk = update_index_lists(ads_index, bulk_index, ec_index, f_ind, b_ind, line,
                                                                                                                           reaction_types_ads,reaction_types_ec,reaction_types_bulk,True, len(species_list[0]))
        elif "-" in line:
            f,b=line[2:].split("-")
            f_ind=index_list(f, species_list)
            b_ind=index_list(b, species_list)
            ads_index, bulk_index, ec_index, reaction_types_ads,reaction_types_ec,reaction_types_bulk = update_index_lists(ads_index, bulk_index, ec_index, f_ind, b_ind, line,
                                                                                                                           reaction_types_ads,reaction_types_ec,reaction_types_bulk,False, len(species_list[0]))
        else:
            print("Mechanism error: cannot index mechanism")
   
    return species_list, [ads_index, bulk_index, ec_index], [reaction_types_ads, reaction_types_bulk, reaction_types_ec]

####################### Functions for generating potential programs

def V_potencial(Ei,Ef,v,amplitude,frequency,nt,f):
    # calculate the potentail input signal for a CV or ACV
    # nt=int(2*nt*abs((f*(Ei)-f*(Ef))))
    En=(Ei+Ef)/2
    tp=2*(Ei-Ef)/v
    a=Ei-En
    ts=abs(2*(-Ei+Ef)/v)
    t=np.linspace(0,ts,nt)
    Edc=a-2*a/np.pi*np.arccos(np.cos(2*np.pi/tp*t))+En+amplitude*np.sin(2*np.pi*frequency*t)
    return Edc,t

############################## Functions for precalc and simulator

def get_kinetic_constants(k_vector, kinetic_types):
    # A function wich checks the number of constants for reversible or irreversible steps and makes the lists correct (assigns a zero for the back step in irreversible steps)
    for i in range(len(k_vector)):
        if kinetic_types[i]==True and len(k_vector[i])!=2:
            print("Error in number of constants: Reversible step, assigned 1 constants, requres 2")
            sys.exit()
            
        elif kinetic_types[i]==False and len(k_vector[i])==1:
            k_vector[i].append(0)
            
        elif kinetic_types[i]==False and len(k_vector[i])!=1:
            print("Error in number of constants: Irreversible step, assigned 2 constants, requres 1")
            sys.exit()
    
    return k_vector

def iterate_Over_conc(step, c, term):
    for i in step:
        term*=c[i]
    return term

def update_K_Matrix(k_matrix, term_f, term_b, step):
    check=[]
    for i in step:
        if i not in check:
            check.append(i)
            k_matrix[i]+=-term_f
            k_matrix[i]+=term_b
    return k_matrix
        
def calc_kinetics(reac_type,c,index,kinetic_const):
    # A function that given the reaction type 0-ads, 1- bulk 
    # the relevant concentrations, indexes connecting the c and constants, 
    # evaluates the reaction forward and backward kinetic rate
    k_matrix=np.zeros(len(c))

    for i in range(len(index[reac_type])):
        
        step=index[reac_type][i]
        constants=kinetic_const[i]
        
        forward_step=constants[0]
        backward_step=constants[1]
        
        forward_step=iterate_Over_conc(step[0], c, forward_step)
        backward_step=iterate_Over_conc(step[1], c, backward_step)
        
        k_matrix=update_K_Matrix(k_matrix, forward_step, backward_step, step[0])
        k_matrix=update_K_Matrix(k_matrix, -forward_step, -backward_step, step[1])
                
    return k_matrix

def calc_EC_kinetics(reac_type, c, index, kinetic_const, E):
    # A function that given the reaction type 0-ads, 1- bulk 
    # the relevant concentrations, indexes connecting the c and constants, 
    # evaluates the reaction forward and backward electrochemical kinetic rate at the boundary
    k_matrix=np.zeros(len(c))

    for i in range(len(index[reac_type])):
        
        step=index[reac_type][i]
        constants=kinetic_const[i]
        
        forward_step=constants[0](E)
        backward_step=constants[1](E)
        
        forward_step=iterate_Over_conc(step[0], c, forward_step)
        backward_step=iterate_Over_conc(step[1], c, backward_step)
        
        k_matrix=update_K_Matrix(k_matrix, forward_step, backward_step, step[0])
        k_matrix=update_K_Matrix(k_matrix, -forward_step, -backward_step, step[1])
                
    return k_matrix

def calc_current(reac_type, c, index, kinetic_const, E):
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
        forward_step=iterate_Over_conc(step[0], c, forward_step)
        backward_step=iterate_Over_conc(step[1], c, backward_step)
        current+=-forward_step+backward_step
    return current

def find_gama(dx,xmax,nx):
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

def Fornberg_weights(z,x,n,m):
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

def Space_ranges(tmax,f,D,fraction,nx):
    # Given the simulation time, f, the maximum diffusion coefficient, the initial dx
    # and the lenghth of spatial direction
    # evaluates a one dimensional grid to be used in simulation
    # fraction is given as dx/xmax
    xmax=6*np.sqrt(tmax*D)
    dx=fraction*xmax
    gama=find_gama(dx, xmax, nx)
    N=np.arange(nx+2)
    x=dx*(gama**N-1)/(gama-1)
    return x

def calc_main_coef(x,dt,D,nx,B):
    # calculate alfas and a's used in simulation
    # calculated with given spatial direction x
    # the weights are given via the method of finite difference implicit method
    # B is to be implemented
    a1=[]
    a2=[]
    a3=[]
    a4=[]
    
    for i in range(1,nx):
        
        weights=Fornberg_weights(x[i],x[i-1:i+3],4,2)
        
        alfa1=weights[0,2]-(B*x[i]**2)*weights[0,1]
        alfa2=weights[1,2]-(B*x[i]**2)*weights[1,1]
        alfa3=weights[2,2]-(B*x[i]**2)*weights[2,1]
        alfa4=weights[3,2]-(B*x[i]**2)*weights[3,1]
        
        a1.append(-alfa1*D*dt)
        a2.append(-alfa2*D*dt+1)
        a3.append(-alfa3*D*dt)
        a4.append(-alfa4*D*dt)
    
    return np.array([np.array(a1),np.array(a2),np.array(a3),np.array(a4)])

def calc_boundary_condition(x,i,D,nx,B):
    # A function for evaluation of the flux boundary condition, at either boundary
    # i should be 0 or -1, 0 for the electrode, -1 for the bulk limit
    # B is used in case of rotation (to be implemented)
    a1=[]
    a2=[]
    a3=[]
    
    if i==0:
        weights=Fornberg_weights(x[i],x[i:i+3],3,1)
    elif i==-1:
        weights=Fornberg_weights(x[i],x[i-2:],3,1)
    else:
        print("Boundary Error: boundary flux indexed incorrectly")
        
    alfa1=weights[0,1]-(B*x[i]**2)
    alfa2=weights[1,1]-(B*x[i]**2)
    alfa3=weights[2,1]-(B*x[i]**2)
    
    a1.append(-alfa1*D)
    a2.append(-alfa2*D)
    a3.append(-alfa3*D)
    
    return np.array([np.array(a1),np.array(a2),np.array(a3)])

def Butler_volmer_kinetics(alpha, k0, E0, f):
    # A function for evaluating the butler-volmer kinetics 
    # it transforms the given constants into function to be evaluated during simulation
    return [lambda E: k0*np.exp(-alpha*f*(E-E0)), lambda E:k0*np.exp((1-alpha)*f*(E-E0))]

def get_EC_kinetic_constants(k_vector, kinetic_types, f):
    # A function for getting BV kinetics at the boundary condition
    # in case of irreversible kinetics the function is a zero function
    for i in range(len(k_vector)):
        k_vector[i]=Butler_volmer_kinetics(k_vector[i][0], k_vector[i][1], k_vector[i][2], f)
        if kinetic_types[i]==False:
            k_vector[i][1]=lambda E: 0
    return k_vector

def time_step(c, a, cp, nx, dt, n1, n, bound1, bound2, pnom, constants, index, F, delta):
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

    bound_kinetics=calc_kinetics(0, np.append(theta, c[0,:]), index, constants[0])+calc_EC_kinetics(2,np.append(theta, c[0,:]), index, constants[2], p)

    f[:n1]=theta-thetap-dt*bound_kinetics[:n1]
    
    f[n1:n1+n]=np.sum(bound1[:,0,:]*c[0:3,:],axis=0)-bound_kinetics[n1:]

    if n!=0:
        for xx in range(1,nx):
            f[n1+n*xx:n1+n*xx+n]=(np.sum(a[:,xx-1,:]*c[xx-1:xx+3,:],axis=0)
                                  -dt*calc_kinetics(1, c[xx,:], index, constants[1])-cp[xx,:])
        f[-2*n:-n]=(c[-2,:])-bound2
        f[-n:]=(c[-2,:]-c[-1,:])
    else:
        pass
        
    ga=F*A*calc_current(2, np.append(theta, c[0,:]), index, constants[2], p)
    
    f9=(1+Ru*Cdl/dt)*gc-Cdl*delta-Ru*Cdl*(gcp)/dt
    f10=pnom-p-Ru*ga-Ru*gc
    f=np.append(f,np.array([f9,f10]))
    return f

def simulator_Main_loop(Mechanism, Constants, Spatial_info, Time, Species_information, Potential_program):
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
    
    F=96485
    R=8.314

    spec, index, types=Parse_mechanism(Mechanism)
    n=len(spec[1])
    n1=len(spec[0])
    
    ads_const, bulk_const, EC_const, cell_const, Diffusion_const = Constants
    
    T,Ru,Cdl,A=cell_const
    f=F/R/T
    
    ads_const=get_kinetic_constants(ads_const, types[0])
    bulk_const=get_kinetic_constants(bulk_const, types[1])
    EC_const=get_EC_kinetic_constants(EC_const, types[2], f)
    
    
    dt=np.average(np.diff(Time))
    
    if len(Species_information[1])>0:
        x=Space_ranges(Time[-1], f, max(Diffusion_const), Spatial_info[0], Spatial_info[1])
    else:
        x=Space_ranges(Time[-1], f, 1, Spatial_info[0], Spatial_info[1])
    
    a=calc_main_coef(x, dt, Diffusion_const, len(x)-2, 0.0)
    
    theta=np.array(Species_information[0])
    
    c=np.zeros((len(x),n))
    for i in range(len(spec[1])):
        c[:,i]=Species_information[1][i](x)
    bound2=c[-1,:] 
    bound1=calc_boundary_condition(x, 0, Diffusion_const, 3, 0)
    
    
    c=c.reshape((1,n*len(x)))[0,:]
    c=np.append(theta,c)
    
    current=[0]
    cap_cur=[0]
    c=np.append(c,np.array([0.25,0]))
    ps=[Potential_program[0]]
    delta_E=np.diff(Potential_program)/dt
    
    constants=[ads_const, bulk_const, EC_const, cell_const]

    for tt in range(1,len(Time)):
        cp=c
        cp[-2]=Potential_program[tt]
        res=sciop.root(time_step, cp, args=(
            a,cp,len(x)-2, dt , n1, n, 
            bound1, bound2, Potential_program[tt], 
            constants, index, F, delta_E[tt-1]),tol=10**-20)

        c=res.x
        current.append(F*A*calc_current(2, c[:n1+n], index, EC_const, c[-2]))
        cap_cur.append(c[-1])
        ps.append(c[-2])
    
    return np.array(ps), np.array(current)+np.array(cap_cur), Time

###################### Functions for initial conditions

def uniform(c):
    # Given a concentration value of c give a uniform function on x
    return lambda x: c


##################### Functions for treating and ploting data
def FFT_analysis(a,f,N,w):
    ### input parameters:
    #a-ndarray of data from voltammetric experimet in form E,I,t
    #f-base frequency of ACV potential
    #N-number of harmonics to be processed
    #w-array or list of window function parameters
    
    ### output parameters:
    #sp-ndarray of FT values of the input current
    #freq-ndarray of frequencies; to be used in a window function or as x-axis
    #I_harmonics-list of ndarrays of analytical signals of the harmonic currents; already procesed
    
    def rectangular(f,w0,w):
        return np.where(abs(f-w0)<=w,1,0)
    
    I_harmonics=[]
    dt=np.average(np.diff(a[:,2]))
    freq=np.fft.fftfreq(a[:,2].shape[-1],d=dt)
    sp=np.fft.fft(a[:,1])
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

def Harmonic_plots(I_harmonics,x_axis,w=0,dt=0,label="",Title=""):
    ###input parameters:
    #I_harmonics-list of ndarrays of analytical harmonic currents to be plotted
    #x_axis-a list or ndarray to plot currents against
    #can be either time or experimental potential
    #w-width of rectangular window function to be used when plotting against potential, dc potential found in FT of experimental potential
    #dt-time step used in experiment; used when finding potential to plot against
    #label-label data in the plot
    #Title-Adding to the title of the plot, format is: "nth harmonic+Title"
    
    
    def rectangular(f,w0,w):
        return np.where(abs(f-w0)<=w,1,0)
    
    if w==0:
        for i in range(len(I_harmonics)):
            plt.figure(str(i)+"_harmonic")
            plt.title(str(i)+"th harmonic "+Title)
            plt.plot(x_axis,I_harmonics[i],label=label)
            plt.xlabel("time [s]")
            plt.ylabel("I [A]")
            # plt.savefig(str(i)+"_harmonic")
    
    elif w<0:
        print("Error in ploting function, w < 0")
    
    elif dt<=0:
        print("Error in ploting function, dt <= 0")
    
    elif w>0 and dt>0:
        sp=np.fft.fft(x_axis)
        freq=np.fft.fftfreq(x_axis.shape[-1],d=dt)
        E_sim=np.fft.ifft(rectangular(freq,0,w)*sp).real
        for i in range(len(I_harmonics)):
            plt.figure(str(i)+"_harmonic")
            plt.title(str(i)+"th harmonic "+Title)
            plt.plot(E_sim,I_harmonics[i],label=label)
            plt.xlabel("$E_{DC}$ [V]")
            plt.ylabel("I [A]")
            # plt.savefig(str(i)+"_harmonic")

def FT_plot(freq,sp,Title="",label=""):
    plt.figure("Fourier_Transform")
    plt.title("Fourier Transform "+Title)
    plt.plot(freq,np.log10(np.abs(sp)),label=label)
    plt.xlabel("frequency [Hz]")
    plt.ylabel("$log_{10}$($I_{FT}$) [dB]")
    
def Plot_measure(x_axis,y_axis,Title="",label="",x_label=""):
    plt.figure("Base_signal")
    plt.title("Measured Signal "+Title)
    plt.plot(x_axis,y_axis,label=label)
    plt.xlabel(x_label)
    plt.ylabel("I [A]")