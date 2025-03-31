from electrokitty import ElectroKitty
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
import time

def plot_chain(chain):
    n_par = len(chain[0, :])
    fig, axes = plt.subplots(n_par)
    for i in range(n_par):
        axes[i].plot(chain[:,i])

if __name__ == "__main__":
    t1 = time.time()
    mechanism = "E(1): c*=a_inb*"

    kin = [[0.5, 1, 0.]]

    D=[]

    spec_info = [[0, 1.11e-05], []]

    si = [0.001, 20, 10**-5, 0]

    cell_const = [293, 0, 0, 0.283*10**-4]

    iso = [0, 0]

    sim = ElectroKitty(mechanism)

    e, t = sim.V_potential(-0.5, 0.5, 0.1, 0, 0, 400)
    sim.create_simulation(kin, cell_const, D, iso, si, spec_info)
        
    e, i, t = sim.simulate()
    i += np.random.normal(0, 0.01*max(i), len(i))
    sim.set_data(e, i, t)
    sim.simulate()
    sim.Plot_data()
    sim.Plot_simulation()
    sim.sample_parameter_distribution(n_samples = 40000, num_chains=2, multi_processing=True, n_processes=2)
    print(len(sim.chains))

    n_par = len(sim.chains[0][0, :])
    fig, axes = plt.subplots(n_par)
    for j in range(len(sim.chains)):
        chain = sim.chains[j]
        for i in range(n_par):
            axes[i].plot(chain[:,i])

    plt.show()
    #plt.plot(e, i/v)