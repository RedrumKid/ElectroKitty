from electrokitty import ElectroKitty

problem = ElectroKitty()

problem.create_simulation([[0.5, 10, 0]], [293, 0, 0, 10**-4], 2*[10**-9], [], [0.0001, 10, 10**-5, 0], [[], [10, 0]])
problem.V_potential(0.5, -0.5, 0.05, 0, 0, 1000)
problem.simulate()
problem.Plot_simulation()

dict = problem.save_json("justatest")
