import lib
import math
from matplotlib import pyplot as plt

# random seed and other parameters
Xo = 1
# (ii)
m2 = 16381
a2 = 572
c = 0
N = 100000
# To store random numbers
randomNums = [0]*N
randomNums2 = lib.mlcg(Xo, m2, a2, c, randomNums,	N)
randomNums2 = [number / (m2-1) for number in randomNums2]

'''
# (B) solving the integral given below by Monte Carlo
#  x^2+y^2=r^2
 y^2+z^2=r^2
# x	=	+/-z
y	=	+/-sqrt(r^2-z^2)
# V_2(r,r)=int_(-r)^r(2sqrt(r^2-z^2))^2dz=(16)/3r^3
'''


def func(x):
    f = (2*math.sqrt(1-x**2))**2
    return f


a = -1
b = 1
ss2 = lib.montecarlo_ran_nums(randomNums2, N, func, a, b)
print('Volume of Steinmetz solid by MC method:',
      "\n For a=572, m=16381 : ", ss2[0])

'''
Volume of Steinmetz solid by MC method:
 For a=572, m=16381 :  5.333646576997727 == 16/3
[Finished in 0.911s]
'''
