import numpy as np
import itertools

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

def check_type(param):
        if type(param) is float or type(param) is int:
            return True
        else:
            return False

def create_parameter_guess(kin, species_information, cell_const, isotherm, fit_kin=True,
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
    
    tells: list containing information that tells the simulator to update a certain parameter
    the first element is the number of sets of kinetic parameters, other values are used to correctly distrubute the guess in the simulator
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
                if check_type(param):
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
        if check_type(cell_const[1]):
            tells.append([2, n, n+1, 1, False])
            initial_guess.append(cell_const[1])
            n+=1
        else:
            tells.append([2, n, n+2, 1, True])
            initial_guess.append(cell_const[1][0].return_mean())
            initial_guess.append(cell_const[1][0].return_sigma())
            n+=2

    if fit_Cdl:
        if check_type(cell_const[2]):
            tells.append([2, n, n+1, 2, False])
            initial_guess.append(cell_const[2])
            n+=1
        else:
            tells.append([2, n, n+2, 2, True])
            initial_guess.append(cell_const[2][0].return_mean())
            initial_guess.append(cell_const[2][0].return_sigma())
            n+=2
        
    if fit_A:
        if check_type(cell_const[3]):
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
            if check_type(species_information[0][ind]):
                scratch_list = [1, n, n+1, ind, False]
                gama_guess.append(species_information[0][ind])
                tells.append(scratch_list)
                n += 1
            else:
                scratch_list = [1, n, n+2, ind, True]
                gama_guess.append(species_information[0][ind][0].return_mean())
                gama_guess.append(species_information[0][ind][0].return_sigma())
                tells.append(scratch_list)
                n += 2
        
        initial_guess+=gama_guess
    
    if fit_iso:
        
        iso_guess = []

        for ind in range(len(isotherm)):
            if check_type(isotherm[ind]):
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

def unpack_fit_params(guess, tells, gamma_position, kin, species_information, cell_const, isotherm):
    """
    Function takes the guess, tells and gammma_position to reconstruct the lists for the simulator
    """
    guess=guess.tolist()
    kinetics=kin
    cell_params=cell_const
    spec_info=species_information
    iso = isotherm
    
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

def create_lower_upper_bounds(guess, tells, potential):
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

        if info[0] != 0:
            if info[-1]:
                lower_bound.append(0)
                lower_bound.append(10**-4)
                upper_bound.append(100*guess[info[1]])
                if guess[info[1]+1] == 0:
                    upper_bound.append(1)
                else:
                    upper_bound.append(100*guess[info[1]+1])
            else:
                lower_bound.append(0)
                upper_bound.append(100*guess[info[1]])

    lower_bound.append(0)
    upper_bound.append(1)

    return lower_bound, upper_bound


kin = [[1, 1, [gaussian_distribution(3, 100)]],
       [2],
       [3, 3],[4,4,4]]

spec_info = [[5,5,5,5], [10,10, 3 , 7]]

cell_c = [0, 10, 11, 12]

iso = [6,6,6]

""" a, b, c = create_parameter_guess(kin, spec_info, cell_c, iso, fit_A=True, 
                                 fit_Cdl=True, fit_gamamax=True, fit_iso=True, fit_kin=False, fit_Ru=True, N_disp=20)

print(a)
print()

print(b)
print()

print(c)
print()

k, cc, s, i = unpack_fit_params(a, b, c, kin, spec_info, cell_c, iso)

print(k)
print()
print(cc)
print()
print(s)
print()
print(i)
print()

lb, ub = create_lower_upper_bounds(a, b, [-0.5, 0.5])

print(lb)
print()
print(ub) """

print(kin.copy())