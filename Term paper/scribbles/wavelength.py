import sympy as sym
import numpy as np
import math
import cmath
import matplotlib.pyplot as plt


# sym.init_printing(use_unicode=False, wrap_line=False)
z = sym.Symbol('z', real=True)
n1 = 2.5
n2 = 1.5
ep0 = 8.854187e-12
m = 1
th = 0
lmdd = np.linspace(550e-9, 650e-9, 10)
Hmin = []
betaa = []
for iijj in range(0, len(lmdd)):
    lmd = lmdd[iijj]
    u = 1.2566e-6
    c = 3e8
    w = 2*math.pi*c/lmd
    d1 = 70e-9
    d2 = 150e-9
    d = 220e-9
    L = 440e-9
    # epsilon with perturbation
    ep1 = 1/2*ep0*(n1**2+n2**2)+1/2*ep0*(n2**2-n1**2)
    ep2 = 1/2*ep0*(n1**2+n2**2)-1/2*ep0*(n2**2-n1**2)
    #  kappa
    k = (1j*(1-math.cos(m*math.pi))*math.sqrt(2)*(n2**2-n1**2) *
         math.cos(2*th))/(2*m*lmd*math.cos(th)*math.sqrt(n2**2+n1**2))
    #  delta beta
    db = 2*(math.sqrt((n1**2+n2**2)/2)*(w/c))*math.cos(th)-m*(2*math.pi/d)
    s = cmath.sqrt(np.conj(k)*k-(db/2)**2)
    #  amplitude
    A = sym.exp(1j*db/2*z)*(s*sym.cosh(s*(L-z))+1j*db/2*sym.sinh(s*(L-z))) / \
        (s*sym.cosh(s*L)+1j*db/2*sym.sinh(s*L))

    # f = open('hcap2.txt', 'w')
    # ff = open('beta2.txt', 'w')

    Ns = 200
    Harray = []
    bb = np.linspace(2*math.pi/lmd, 4*math.pi/lmd, 20)
    hmin = 10e22
    for kk in range(0, len(bb)):
        b = bb[kk]
        # f = A*sym.exp(1j*b*z)
        chi = A*sym.exp(1j*b*z)
        chi_chi = np.conj(chi)*chi
        # probability distribution function
        pr = np.conj(chi)*chi/sym.integrate((chi_chi), (z, 0, L))
        # print(pr.subs({z: 1.51347572901559}).evalf())
        Ev = 0
        zz = 0
        step = 0.1e-9
        for ii in range(1, Ns):
            prv = pr.subs({z: zz}).evalf()
            if (0 <= zz and zz <= d1) or (d <= zz and zz <= d1+d):
                # print("zz1", zz)
                H = (sym.diff(chi, z, 2)+w**2*u*(ep1)*chi)
                E = (1/chi)*H*chi
                Ev += E.subs({z: zz}).evalf()
            else:
                # print("zz2", zz)
                H = (sym.diff(chi, z, 2)+w**2*u*(ep2)*chi)
                E = (1/chi)*H*chi
                Ev += E.subs({z: zz}).evalf()
            zzn = zz + np.random.rand()*step
            prv1 = pr.subs({z: zzn}).evalf()
            R = np.abs(prv1/prv)
            if R < 1:
                neta = np.random.rand()
                if neta < R:
                    zz = zzn
            else:
                zz = zzn
            if zz >= L:
                zz = 0
        Hcap = Ev/Ns
        # f.write(str(np.abs(Hcap)))
        # f.write('\n')
        # ff.write(str(b))
        # ff.write('\n')
        print(np.abs(Hcap))
        Harray.append(np.abs(Hcap))
        if np.abs(Hcap) < hmin:
            hmin = np.abs(Hcap)
            beta = b
    # f.close()
    # ff.close()

    print(Harray)
    Hmin.append(hmin)
    betaa.append(beta)
    # fig, ax = plt.subplots(figsize=(12, 6))
    plt.plot(bb, Harray, label=lmd)
    plt.legend()
    # plt.savefig('book_read.png')
figg, axx = plt.subplots(figsize=(12, 6))
axx.plot(lmdd, betaa, label='vs')
print("lmadd", lmdd)
print("betaa", betaa)
axx.legend()
plt.show()
