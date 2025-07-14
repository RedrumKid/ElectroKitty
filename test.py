from electrokitty import ElectroKitty
import numpy as np

""" mechanism = "E(1): a* < b* \n C: b* = c* \n C: a > b + b \n C: c>l"

kin = [[0.5, 10, 0]]

ic = [[10**-5, 0], []]

D = []

iso = []

si = [0.001, 10, 10**-5, 0]

cc = [293, 0, 0, 10**-4]

sim = ElectroKitty(mechanism)

sim.C_potential(0, 100, 1000)

sim.create_simulation(kin, cc, D, iso, si, ic)

spec, index, r_t, r_ind, es = sim.mechanism_list
print(r_ind)
iso = [
    [[1,1,1], [2,2,2]],
    [[3,3,3], [4,4,4]],
    [[5,5,5], [6,6,6]],
    [[7,7,7], [8,8,8]]
]

def update_iso(iso, spec):
    for step in iso:
        for li in step:
            li += len(spec[1])*[0]
    
    return iso

print(update_iso(iso, spec))

def iterate_Over_conc(step, term, isotherm):
        
        for i in step:
            term+=isotherm[i]
        return term


def calc_kinetics(reac_type, index, r_ind, isotherm):
    
        for i in range(len(index[reac_type])):
            
            step=index[reac_type][i]
            iso = isotherm[r_ind[reac_type][i]]
            forward_step=0
            backward_step=0
            print(iso)
            forward_step=iterate_Over_conc(step[0], forward_step, iso[0])
            backward_step=iterate_Over_conc(step[1], backward_step, iso[1])

            print(forward_step)
            print(backward_step)

        return forward_step, backward_step

calc_kinetics(2, index, r_ind, iso) """


def check(param):
    if type(param) is not list and type(param) is not tuple and type(param) is not np.array:
        return True
    else:
        return False

a = np.array([5.])

print(check(a))

print(a)
