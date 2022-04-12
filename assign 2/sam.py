import lib
import numpy as np
import matplotlib.pyplot as plt

opn = open('B2.txt', 'r')
lsplit = opn.readline().split()
b = []
for val in lsplit:
    b.append(float(val))
rows = 6
with open('A2.txt') as f:
    a = []
    for i in range(0, rows):
        a.append(list(map(float, f.readline().split())))

with open('Bi.txt') as f:
    bi = []
    for i in range(0, rows):
        bi.append(list(map(float, f.readline().split())))

# solution using jacobi
aj = np.array(a)
bj = np.array(b)
bij = np.array(bi)
ITERATION_LIMIT = 100000
j = lib.jacobimat(aj, bj, ITERATION_LIMIT)

print('Solution by jacobi :', j[0])
# print("Error:", j[1])

# Solution using LU
# x or solution of given linear equations
x = lib.forback(a, b)[1]
print('\n Solution by LU :', x)


tol = 1e-4
# inverse using gaussjordan
invr = lib.inverse(a)
print('\n Inverse :', invr)

# inverse using gauss_seidel
ings = np.zeros((len(a), len(a)))  # Pre-allocate matrix
for i in range(len(bi)):
    zz = lib.gauss_seidel(aj, bij[i], tol)
    for j in range(1, len(bi)):
        ings[:, i] = zz[0]
print('\n Inverse using Gauss Seidel method :\n', ings)

# inverse using Conjugate Gradient
incg = np.zeros((len(a), len(a)))
x0 = np.ones(len(a))
for i in range(len(bi)):
    cg = lib.conjGrad(aj, bij[i], x0, tol)
    for j in range(1, len(bi)):
        incg[:, i] = cg[0]
print('\n Inverse using Conjugate gradient method :\n', incg)

# inverse using jacobi
injc = np.zeros((len(a), len(a)))
x0 = np.ones(len(a))
for i in range(len(bi)):
    jc = lib.conjGrad(aj, bij[i], x0, tol)
    for j in range(1, len(bi)):
        injc[:, i] = jc[0]
print('\n Inverse using Jacobi method :\n', injc)


# define subplot grid
fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(15, 12))
plt.subplots_adjust(hspace=0.5)
fig.suptitle("Residue versus iteration steps", fontsize=18, y=0.95)
# gauss_seidel
axs[1].plot(zz[2], zz[1])
axs[1].set_title("Gauss Seidel")
axs[1].set_ylabel('Residue')
axs[1].set_xlabel('No. of iteration')
# Conjugate gradient
axs[0].plot(cg[2], cg[1])
axs[0].set_title("Conjugate Gradient")
axs[0].set_ylabel('Residue')
axs[0].set_xlabel('No. of iteration')
# Jacobi
axs[2].plot(jc[2], jc[1])
axs[2].set_title("Jacobi")
axs[2].set_ylabel('Residue')
axs[2].set_xlabel('No. of iteration')
plt.show()


'''
System:
2.0*x1 + -3.0*x2 + 0.0*x3 + 0.0*x4 + 0.0*x5 + 0.0*x6 = -1.66666
-1.0*x1 + 4.0*x2 + -1.0*x3 + 0.0*x4 + -1.0*x5 + 0.0*x6 = 0.66666
0.0*x1 + -1.0*x2 + 4.0*x3 + 0.0*x4 + 0.0*x5 + -1.0*x6 = 3.0
0.0*x1 + 0.0*x2 + 0.0*x3 + 2.0*x4 + -3.0*x5 + 0.0*x6 = -1.33333
0.0*x1 + -1.0*x2 + 0.0*x3 + -1.0*x4 + 4.0*x5 + -1.0*x6 = -0.333333
0.0*x1 + 0.0*x2 + -1.0*x3 + 0.0*x4 + -1.0*x5 + 4.0*x6 = 1.666666

Solution by jacobi : [-3.33332182e-01  3.33331879e-01  9.99999576e-01
 -6.66664818e-01  1.21099390e-07  6.66666424e-01]
 Solution by LU :   [-0.33333218181818175, 0.33333187878787884,
 0.9999995757575757, -0.666664818181818, 1.2121212131308088e-07, 0.6666664242424243]

 Inverse using Gauss Seidel method :
 [[0.93490644 0.86998776 0.2596891  0.20766252 0.41541074 0.16876865]
 [0.28997254 0.58002315 0.17313734 0.13847025 0.27697875 0.11252621]
 [0.08655462 0.17313734 0.3203381  0.05625622 0.11252621 0.10821506]
 [0.20766252 0.41546813 0.16878931 0.93495881 0.86998776 0.2596891 ]
 [0.13847025 0.27700437 0.11253543 0.28999592 0.58002315 0.17313734]
 [0.05625622 0.11253543 0.10821838 0.08656303 0.17313734 0.3203381 ]]

 Inverse using Conjugate gradient method :
 [[0.93508434 0.87017068 0.25974643 0.20776419 0.41555179 0.16887971]
 [0.29006418 0.58010308 0.17316788 0.13852489 0.27705023 0.11257167]
 [0.08658568 0.17317276 0.32036005 0.0562844  0.11255597 0.10822718]
 [0.20776419 0.41555179 0.16887971 0.93508434 0.87017068 0.25974643]
 [0.13852489 0.27705023 0.11257167 0.29006418 0.58010308 0.17316788]
 [0.0562844  0.11255597 0.10822718 0.08658568 0.17317276 0.32036005]]

 Inverse using Jacobi method :
 [[0.93508434 0.87017068 0.25974643 0.20776419 0.41555179 0.16887971]
 [0.29006418 0.58010308 0.17316788 0.13852489 0.27705023 0.11257167]
 [0.08658568 0.17317276 0.32036005 0.0562844  0.11255597 0.10822718]
 [0.20776419 0.41555179 0.16887971 0.93508434 0.87017068 0.25974643]
 [0.13852489 0.27705023 0.11257167 0.29006418 0.58010308 0.17316788]
 [0.0562844  0.11255597 0.10822718 0.08658568 0.17317276 0.32036005]]
[Finished in 3.586s]

'''
