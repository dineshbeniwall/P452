import lib

opn = open('B1.txt', 'r')
lsplit = opn.readline().split()
b = []
for val in lsplit:
    b.append(float(val))
rows = 6
with open('A1.txt') as f:
    a = []
    for i in range(0, rows):
        a.append(list(map(float, f.readline().split())))

# Gausss Jordan
lib.gaussjordan(a, b)
print('Solution by Gauss jordan : \n matrix', b)


# LU decomposion
# upper triangle matrix
ut = lib.ludecom(a)[1]
# lower triangle matrix
lt = lib.ludecom(a)[0]
# AX=B      LUX=B   LZ=B   UX=Z
# Z
z = lib.forback(a, b)[0]
# x or solution of given linear equations
x = lib.forback(a, b)[1]
print('Solution by LU : \n matrix ', x)

'''

Solution by Gauss jordan :
 matrix [-1.7618170439978567, 0.8962280338740136, 4.051931404116157, -1.6171308025395428, 2.041913538501914, 0.15183248715593495]
Solution by LU :
 matrix  [-1.7618170439978567, 0.8962280338740136, 4.051931404116157, -1.6171308025395428, 2.041913538501914, 0.15183248715593495]
[Finished in 0.717s]

'''
