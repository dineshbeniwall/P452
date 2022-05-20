import math
import lib
from matplotlib import pyplot as plt
import numpy as np

# random seed and other parameters
Xo = 1
# (i)
m1 = 16381
a1 = 572
c = 0
N = 200
M = 500
# To store random numbers
randomNums = [0]*N
# Function to generate random numbers
# multiplicative linear congruential generator
randomNums1 = lib.mlcg(Xo, m1, a1, c, randomNums,	N)
# scale them to [0,1]
randomNums1 = [number / (m1-1) for number in randomNums1]
randomNums1 = np.array(randomNums1)
# scale them to [0,2*pi]
randomNums1 = [number*2*math.pi for number in randomNums1]
# print(randomNums1[10])


# random walk rnwanalysis
# a = lib.randomwalk(100)
aa = lib.mod_rnwanalysis(N, M, Xo, m1, a1, c, randomNums)[0]
rms = lib.mod_rnwanalysis(N, M, Xo, m1, a1, c, randomNums)[4]
print("RMS : ", rms, "  Sqrt(N) : ", math.sqrt(N))
# plotting first 5 wlaks
aa = np.array(aa, dtype=object)
for i in range(5):
    plt.plot(aa[i-1][0], aa[i-1][1])
plt.title('First 5 walks of 200 steps out of 500')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.show()


# # R_rms(N) = sqrt(N)
# # output of mod_rnwanalysis return aa, Dx, Dy, avgr, rms
NN = []
RMS = []
N = 50

for i in range(10):
    randomNums = [0]*N
    rms = lib.mod_rnwanalysis(N, M, Xo, m1, a1, c, randomNums)[4]
    # a = lib.rnwanalysis(500, 10)[0]
    NN.append(math.sqrt(N))
    RMS.append(rms)
    # print("rms :", rms, "N", math.sqrt(N))
    N += 50

plt.plot(RMS, NN)
plt.title('R_rms(N) vs Sqrt(N)')
plt.xlabel('$sqrt(N)$')
plt.ylabel('$R_rms(N)$')
plt.show()


'''
For N=500 and 500 walks
RMS :  15.363270250839534   Sqrt(N) :  14.142135623730951
[Finished in 47.881s]

'''
