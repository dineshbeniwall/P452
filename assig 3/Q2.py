import numpy as np
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
'''
Use shooting method with RK4 integrator to solve Schr¨odinger
equation for infinite square well potential, of unit width,
for the two lowest energy state.
'''


def Schroed(y, r, V, E):
    psi, phi = y
    dphidx = [phi, (V-E)*psi]
    return np.array(dphidx)


def rk4(f, psi0, x, V, E):
    n = len(x)
    psi = np.array([psi0]*n)
    for i in range(n-1):
        h = x[i+1]-x[i]
        k1 = h*f(psi[i], x[i], V[i], E)
        k2 = h*f(psi[i]+0.5*k1, x[i]+0.5*h, V[i], E)
        k3 = h*f(psi[i]+0.5*k2, x[i]+0.5*h, V[i], E)
        k4 = h*f(psi[i]+k3, x[i+1], V[i], E)
        psi[i+1] = psi[i]+(k1+2.0*(k2+k3)+k4)/6.0
    return psi


def findZeros(rightbound_vals):
    return np.where(np.diff(np.signbit(rightbound_vals)))[0]


def normalize(output_wavefunc):
    normal = max(output_wavefunc)
    return output_wavefunc*(1/(normal))


def countNodes(wavefunc):
    maxarray = argrelextrema(wavefunc, np.greater)[0]
    minarray = argrelextrema(wavefunc, np.less)[0]
    nodecounter = len(maxarray)+len(minarray)
    return nodecounter


def RefineEnergy(Ebot, Etop, Nodes, psi0, x, V):
    tolerance = 1e-12
    ET = Etop
    EB = Ebot
    psi = [1]
    while (abs(EB-ET) > tolerance or abs(psi[-1]) > 1e-3):
        initE = (ET+EB)/2.0
        psi = rk4(Schroed, psi0, x, V, initE)[:, 0]
        nodes_ist = len(findZeros(psi))-1
        if nodes_ist > Nodes+1:
            ET = initE
            continue
        if nodes_ist < Nodes-1:
            EB = initE
            continue
        if (nodes_ist % 2 == 0):
            if((psi[len(psi)-1] <= 0.0)):
                ET = initE
            else:
                EB = initE
        elif nodes_ist > 0:
            if((psi[len(psi)-1] <= 0.0)):
                EB = initE
            else:
                ET = initE
        elif nodes_ist < 0:
            EB = initE
    return EB, ET


def ShootingInfinitePotentialWell(E_interval, nodes):
    psi_0 = 0.0
    phi_0 = 1.0
    psi_init = np.array([psi_0, phi_0])
    h_mesh = 1.0/100.0
    x_arr_ipw = np.arange(0.0, 1.0+h_mesh, h_mesh)
    V_ipw = np.zeros(len(x_arr_ipw))
    EBref, ETref = RefineEnergy(E_interval[0], E_interval[1], nodes, psi_init, x_arr_ipw, V_ipw)
    psi = rk4(Schroed, psi_init, x_arr_ipw, V_ipw, EBref)[:, 0]
    return EBref, normalize(psi), x_arr_ipw


# Start
E_ipw = [1.0, 500.0]
nodes_arr = [0, 1]
L = 0.0
N = 1.0
figipw = plt.figure()
Energy, psi_ipw, x_ipw = ShootingInfinitePotentialWell(E_ipw, 0)
print('Energy (Ground state): %f ' % Energy)
plt.plot(x_ipw, psi_ipw)
plt.title('Eigenstate : % s' % (0, ))
plt.xlabel('x')
plt.ylabel(r'$\psi$')
plt.show()

Energy, psi_ipw, x_ipw = ShootingInfinitePotentialWell(E_ipw, 2)
print('Energy (First excited state): %f ' % Energy)
plt.plot(x_ipw, psi_ipw)
plt.title('Eigenstate : % s' % (1, ))
plt.xlabel('x')
plt.ylabel(r'$\psi$')
plt.show()

'''
Energy (Ground state): 9.869605
Energy (First excited state): 39.478428
[Finished in 9.748s]
'''
