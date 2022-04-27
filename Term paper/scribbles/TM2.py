import sympy as sym
import numpy as np
import math
import cmath
from scipy.integrate import quad
import lib

# sym.init_printing(use_unicode=False, wrap_line=False)
z = sym.Symbol('z', real=True)
n1 = 2.5
n2 = 1.5
ep0 = 8.854187e-12
m = 1
th = 0
lmd = 600e-9
u = 1.2566e-6
c = 3e8
w = 2*math.pi*lmd/c
d1 = 70e-9
d2 = 150e-9
d = 220e-9
L = 220e-9
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
A = sym.exp(1j*db/2*z)*(s*sym.cosh(s*(d-z))+1j*db/2*sym.sinh(s*(d-z))) / \
    (s*sym.cosh(s*d)+1j*db/2*sym.sinh(s*d))
# print(A.subs({z: 22e-90}).evalf())
b = 1
# f = A*sym.exp(1j*b*z)
chi = A*sym.exp(1j*b*z)
chi_chi = np.conj(chi)*chi

# print(np.abs(sym.integrate(chi_chi, (z, 0, L))))
# print(chi_chi.subs({z: 200e-9}).evalf())
# # probability distribution function
pr = np.abs(chi_chi/sym.integrate((chi_chi), (z, 0, L)))
print(min(np.linspace(-1, 1, 1)))

'''
def f(x):
    return np.abs(chi_chi.subs({z: x}).evalf())


inter = lib.montecarlo(0, 220e-9, 100, f)
print(inter)
'''
'''
Ns = 50
Harray = []
carray = []
ccarray = []
zz = np.linspace(0, 220e-9, 100)
for kk in range(0, len(zz)):
    c = np.abs(chi.subs({z: zz[kk]}).evalf())
    cc = np.abs(chi_chi.subs({z: zz[kk]}).evalf())
    carray.append(c)
    ccarray.append(cc)
# plt.plot(zz, carray)
plt.plot(zz, ccarray)
plt.show()
'''
'''
for kk in range(0, len(bb)):
    b = bb[kk]
    # f = A*sym.exp(1j*b*z)
    chi = A*sym.exp(1j*b*z)
    chi_chi = np.conj(chi)*chi
    # probability distribution function
    pr = np.abs(np.conj(chi)*chi)/sym.integrate((chi_chi), (z, 0, L))
    # print(pr.subs({z: 1.51347572901559}).evalf())
    Ev = 0
    zz = 0.001e-9
    step = 0.000001
    for ii in range(1, Ns):
        prv = pr.subs({z: zz}).evalf()
        if 0 <= zz and zz <= d1:
            print("zz1", zz)
            H = (sym.diff(chi, z, 2)+w**2*u*(ep1)*chi)
            E = (1/chi)*H*chi
            Ev += E.subs({z: zz}).evalf()
        else:
            print("zz2", zz)
            H = (sym.diff(chi, z, 2)+w**2*u*(ep2)*chi)
            E = (1/chi)*H*chi
            Ev += E.subs({z: zz}).evalf()
        zz = zz+prv*step
        prv1 = pr.subs({z: zz}).evalf()
        print("prv 1", prv1)

        if prv1/prv < 1:
            for j in range(1, 1000):
                zz = np.random.rand()
                if zz < prv1/prv:
                    print("zz3", zz)
                    break

    Hcap = Ev/Ns
    # print(np.abs(Hcap))
    Harray.append(np.abs(Hcap))

print(Harray)
plt.plot(bb, Harray)
plt.show()
'''
