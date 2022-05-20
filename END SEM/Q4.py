import math


def LP(n):
    if n == 4:
        r = [0.861136311, 0.339981043, -0.339981043, -0.861136311]
        w = [0.347854845, 0.652145154, 0.652145154, 0.347854845]
    elif n == 5:
        r = [0.906179845, 0.538469310, 0, -0.538469310, -0.906179845]
        w = [0.236926885, 0.478628670, 0.568888889, 0.478628670, 0.236926885]
    elif n == 6:
        r = [0.932469514, 0.661209386, 0.238619186, -0.238619186, -0.661209386, -0.932469514]
        w = [0.171324492, 0.360761573, 0.467913934, 0.467913934, 0.360761573, 0.171324492]
    return (r, w)


def int_GausQuad(f, n):
    int = 0
    (r, w) = LP(n)
    for i in range(n):
        int += w[i-1] * f(r[i-1])
    return int


def fun(x):
    return 1/math.sqrt(1+x**2)


degree = [4, 5, 6]
for n in degree:
    v = int_GausQuad(fun, n)
    print("Potential for degree ", n, ' : ', v)

'''
Potential for degree  4  :  1.7620541789046658
Potential for degree  5  :  1.7628552954010728
Potential for degree  6  :  1.762730048499759
[Finished in 0.422s]
'''
