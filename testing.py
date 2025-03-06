import numpy as np
import itertools

def f(x):
    y = 1
    for val in x:
        if val>0.5:
            y*=1
        else:
            y*=0
    return y

N = 2

y = 1/N*np.ones(N)
x = np.linspace(0, 1, N)

p = lambda x: 1/N

l = [x]
l1 = [y]

suma = 0
mult = 1

suma = 0

""" for combinations in zip(itertools.product(*l), itertools.product(*l1)):
    product = 1
    for value in combinations[1]:
        product *= value
    suma += product*f(combinations[0])

print(suma) """

arr = [[0,1], [1,2]]
print(arr[[0,1]])
#print(np.meshgrid(*l))

