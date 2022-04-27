import sympy as sym
import numpy as np
import math


sym.init_printing(use_unicode=False, wrap_line=False)
n1 = sym.Symbol('n1', real=True)
n2 = sym.Symbol('n2', real=True)
ep0 = sym.Symbol('ep0', real=True)
m = sym.Symbol('m', real=True)
th = sym.Symbol('th', real=True)
lmd = sym.Symbol('lmd', real=True)
z = sym.Symbol('z', real=True)
w = sym.Symbol('w', real=True)
u = sym.Symbol('u', real=True)
c = sym.Symbol('c', real=True)
d = sym.Symbol('d', real=True)
b = sym.Symbol('b')
# epsilon with perturbation
ep = 1/2*ep0*(n1**2+n2**2)+1/2*ep0*(n2**2-n1**2)
#  kappa
k = (1j*(1-sym.cos(m*sym.pi))*sym.sqrt(2)*(n2**2-n1**2) *
     sym.cos(2*th))/(2*m*lmd*sym.cos(th)*sym.sqrt(n2**2+n1**2))
#  delta beta
db = 2*(sym.sqrt((n1**2+n2**2)/2)*(w/c))*sym.cos(th)-m*(2*sym.pi/d)
s = sym.sqrt(sym.conjugate(k)*k-(db/2)**2)
#  amplitude
A = sym.exp(1j*db/2*z)*(s*sym.cosh(s*(d-z))+1j*db/2*sym.sinh(s*(d-z))) /\
    (s*sym.cosh(s*d)+1j*db/2*sym.sinh(s*d))
#  f =  A*exp(1i*b*z);
f = A*sym.exp(1j*b*z)
#H = (sym.diff(f, z, 2)+w**2*u*(ep)*f)
ww = sym.conjugate(f)*f
print(np.linspace(-1, 1, 1))
# wHw_val = wHw.subs({th: 0, m: 1, n1: 1, n2: 2, lmd: 1, ep0: 1, w: 1, u: 1, c: 1, d: 4}).evalf()

# Int = sym.integrate((wHw), (z, 0, 4))
# print(Int)
# print(A.subs({th: 0, m: 1, n1: 1, n2: 2, lmd: 1}).evalf())
