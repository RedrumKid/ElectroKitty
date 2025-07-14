from electrokitty import ElectroKitty
import numpy as np
import matplotlib.pyplot as plt

mechanism = "E(1): a* = b* \n E(1): b* = c* \n C: a = b \n E(1): b = c \n E(1): c + c* = d*"

kin = [[0.5, 1000, 0.4],
       [0.5, 1000, -0.4],
       [10, 100],
       [0.5, 100, 0],
       [0.5, 1000, -0.2]]

ic = [[10**-5, 0, 0, 0], [1, 0, 0]]

D = 3*[10**-9]

iso = [-6, -4, -1, -10]

""" iso = [[[-6, -4, 0, 0], [-6, -4, 0, 0]],
       [[0., -4, -1, 0], [0., -4, -1, 0]],
       [[0, 0, 0, 0], [0, 0, 0, 0]],
       [[0, 0, 0, 0], [0, 0, 0, 0]],
       [[0, 0, -1, -10], [0, 0 , -1, -10]]] """

si = [0.001, 10, 10**-5, 0]

cc = [293, 0, 0, 10**-4]

sim = ElectroKitty(mechanism)

sim.V_potential(0.75, -0.75, 0.1, 0, 0, 1500)

sim.create_simulation(kin, cc, D, iso, si, ic)

l = len(sim.mechanism_list[0][0])
index = sim.mechanism_list[1]
r_ind = sim.mechanism_list[-2]
size = len(r_ind[0]) + len(r_ind[1]) + len(r_ind[2])

iso1 = size*[[l*[0], l*[0]]]

ind = 2

for j in range(len(r_ind[ind])):
    f, b = index[ind][j]
    print(iso1)
    for el in f:
        if el < l:
            print(iso1[r_ind[ind][j]])
            iso1[r_ind[ind][j]][0][el] = iso[el]
            iso1[r_ind[ind][j]][1][el] = iso[el]
    for el in b:
        if el < l:
            
            iso1[r_ind[ind][j]][0][el] = iso[el]
            iso1[r_ind[ind][j]][1][el] = iso[el]

print(iso1)

""" e, i, t = sim.simulate()
sim.Plot_Adsorbed_species()
plt.plot(e, i)
plt.show() """