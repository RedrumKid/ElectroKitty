from electrokitty import ElectroKitty
import numpy as np

problem = ElectroKitty()

a = [[np.float64(0.5), np.float64(0.5), np.float64(0.)]]

problem.create_simulation(a, [293, 0, 0, np.float32(10**-4)], 2*[10**-9], [], [0.0001, 10, 10**-5, 0], [[], [10, 0]])
problem.V_potential(0.5, -0.5, 0.05, 0, 0, 1000)
e, i, t = problem.simulate()
problem.Plot_simulation()
problem.set_data(e, i, t)

# problem.fit_to_data(algorithm = "CMA-ES", tolx = 10**-2, tolf=10**-2)
# problem.sample_parameter_distribution(n_samples = 60)

# problem.FFT_analyze_sim(10, 2, [1, 1, 1, 1])

dict = problem.save_json("justatest")
# print(dict.keys())

p1 = ElectroKitty()
p1.load_from_json("justatest.json")
p1.Plot_simulation()