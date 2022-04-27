import sympy as sym
import numpy as np
import math
import cmath
# for printing
sym.init_printing(use_unicode=False, wrap_line=False)
x = sym.Symbol('x')
y = sym.Symbol('y')
d = sym.Symbol('d', real=True)

print(sym.integrate((sym.exp(x+d+2j)), (x, 0, d)))
in_complx1 = sym.exp(x+d+4j)
print(sym.integrate((in_complx1), (x, 0, 1)))
out_complx1 = np.conj(in_complx1)
print("Output conjugated : ", out_complx1)
b = np.linspace(-1, 1, 100)
print(b[1])
a = 2+3j
d = 2j
if a/d < 1:
    print("yes")
else:
    print("no")

ep0 = 8.854187e-12
m = 1
th = 0
lmd = 650e-9
w = 2e10
u = 1.2566e-6
c = 3e8
