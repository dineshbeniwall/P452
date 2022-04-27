import sympy as sym
import numpy as np
import math
import cmath
import matplotlib.pyplot as plt


def ep1(n1, n2, ep0):
    ep1 = 1/2*ep0*(n1**2+n2**2)+1/2*ep0*(n2**2-n1**2)
    return ep1


def ep2(n1, n2, ep0):
    ep2 = 1/2*ep0*(n1**2+n2**2)-1/2*ep0*(n2**2-n1**2)
    return ep2


def k(m, n1, n2, th, lmd):
    #  kappa
    k = (1j*(1-math.cos(m*math.pi))*math.sqrt(2)*(n2**2-n1**2) *
         math.cos(2*th))/(2*m*lmd*math.cos(th)*math.sqrt(n2**2+n1**2))
    return k


def db(n1, n2, w, c, th, m, d):
    #  delta beta
    db = 2*(math.sqrt((n1**2+n2**2)/2)*(w/c))*math.cos(th)-m*(2*math.pi/d)
    return db


def s(k, db):
    s = cmath.sqrt(np.conj(k)*k-(db/2)**2)
    return s


def A(db, s, d, z):
    A = sym.exp(1j*db/2*z)*(s*sym.cosh(s*(d-z))+1j*db/2*sym.sinh(s*(d-z))) / \
        (s*cmath.cosh(s*d)+1j*db/2*cmath.sinh(s*d))
    return A
