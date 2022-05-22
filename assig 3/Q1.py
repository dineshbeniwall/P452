import lib
import math
import numpy as np


def func(x):
    return np.exp(-x**2)


Xo = 1
m1 = 16381
a1 = 572
c = 0
N = 100000
# To store random numbers
randomNums = [0]*N
randomNums1 = lib.mlcg(Xo, m1, a1, c, randomNums,	N)
randomNums2 = [number / (m1-1) for number in randomNums1]


# without importance sampling
mc = lib.montecarlo_ran_nums(randomNums2, N, func, 0, 1)
print("Montecarlo Integration without importance sampling : ", mc[0])


# with importance sampling
def get_integral_(m, randomNums1, no_of_steps):
    integral = 0
    for i in range(no_of_steps):
        # Xo, m, a, c, randomNums, N
        x = randomNums1[i]  # a random number is generated
        # from an exponential distribution in [0, 1)
        X = -math.log(1-(math.exp(1)-1)/math.exp(1)*x/m)
        integral += (math.exp(1)-1)/math.exp(1)*math.exp(-X**2)/math.exp(-X)
    return integral/no_of_steps


mcim = get_integral_(m1, randomNums1, N)

print("Montecarlo Integration with importance sampling : {}".format(mcim))

'''
Montecarlo Integration without importance sampling :  0.7467599169351836
Montecarlo Integration with importance sampling : 0.7468486922640555
[Finished in 1.473s]
'''
