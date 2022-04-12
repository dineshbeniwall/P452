import numpy as np
import matplotlib.pyplot as plt


# Operator given in the question
# A = 1/2(d(x+u,y)-d(x-u,y)+2d(x,y))+m^2*d(x,y)
def Amat(i, v, m, n):
    x = i % n
    y = i//n
    Ai = 0
    # boundary condition of X
    if(x == 0):
        Ai += 0.5*v[y*n+n-1]
    else:
        Ai += 0.5*v[y*n+x-1]
    if(x == n-1):
        Ai += 0.5*v[y*n]
    else:
        Ai += 0.5*v[y*n+x+1]
    # boundary condition of Y
    if(y == 0):
        Ai += 0.5*v[n*n-n+x]
    else:
        Ai += 0.5*v[y*n-n+x]
    if(y == n-1):
        Ai += 0.5*v[x]
    else:
        Ai += 0.5*v[y*n+n+x]
    return Ai+(m**2)*v[i]


# dot product with A matrix
def multiply(x, m, n, N):
    xAx = 0
    Ax = np.zeros(n)
    for i in range(n):
        Ax[i] = Amat(i, x, m, N)
        xAx += x[i]*Ax[i]
    return xAx, Ax


# Conjugate Gradient
# A has been replaced by X
# m for the input in Amat
# N is for Amat
def conjGrad(x, b, m, N, tol):
    Nitr = 1000
    n = b.shape[0]
    # r = b - A.dot(x) compare it with
    _, Ax = multiply(x, m, n, N)
    r = b-Ax
    # norm_r = r / np.sqrt(np.sum(r**2))
    p = r.copy()
    er = np.array([])
    kk = np.array([])
    for i in range(Nitr):
        # changed
        pAp, Ap = multiply(p, m, n, N)
        alpha = np.dot(p, r)/np.dot(p, Ap)
        x = x + alpha*p
        # changed
        r = b - alpha*Ap
        er = np.append(er, np.sum((r**2)))
        kk = np.append(kk, i)
        if np.sqrt(np.sum((r**2))) < tol:
            # print('Itr:', i)
            break
        else:
            beta = -np.dot(r, Ap)/np.dot(p, Ap)
            p = r + beta*p
    return x, er, kk


# Main functios
# Dimantion
D = 2
N = 30
m = 0.1
tol = 1e-6

# inverse using Conjugate Gradient
incg = np.zeros((N**2, 2))
bij = np.zeros((N**2, 2))
bij[1, 0] = bij[0, 1] = 1
x0 = np.ones(len(incg))
for i in range(2):
    cg = conjGrad(x0, bij[:, i], m, N, tol)
    for j in range(1, len(incg)):
        incg[:, i] = cg[0]
print('\n Inverse using Conjugate gradient method :\n', incg)

# define subplot grid
fig, axs = plt.subplots(figsize=(10, 10))
plt.subplots_adjust(hspace=0.5)
fig.suptitle("Residue versus iteration steps", fontsize=18, y=0.95)
# CG
axs.plot(cg[2], cg[1])
axs.set_title("Conjugate Gradient")
axs.set_ylabel('Residue')
axs.set_xlabel('No. of iteration')
plt.show()
