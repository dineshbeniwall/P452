import sympy as sym
import numpy as np
import math
import cmath


def normalize(b, n1, n2, m, th, lmd, w, c, d, L):
    sym.init_printing(use_unicode=False, wrap_line=False)
    z = sym.Symbol('z', real=True)
    #  kappa
    k = (1j*(1-math.cos(m*math.pi))*math.sqrt(2)*(n2**2-n1**2) *
         math.cos(2*th))/(2*m*lmd*math.cos(th)*math.sqrt(n2**2+n1**2))
    #  delta beta
    db = 2*(math.sqrt((n1**2+n2**2)/2)*(w/c))*math.cos(th)-m*(2*math.pi/d)
    s = cmath.sqrt(np.conj(k)*k-(db/2)**2)
    #  amplitude
    A = sym.exp(1j*db/2*z)*(s*sym.cosh(s*(d-z))+1j*db/2*sym.sinh(s*(d-z))) / \
        (s*cmath.cosh(s*d)+1j*db/2*cmath.sinh(s*d))
    # f = A*sym.exp(1j*b*z)
    chi = A*sym.exp(1j*b*z)
    ww = np.conj(chi)*chi

    # integrate
    Int = sym.integrate((ww), (z, 0, L))
    return Int
