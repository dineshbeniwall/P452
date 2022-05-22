import lib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


x = np.linspace(0, 1, 10)
y = np.linspace(0, 1, 10)

potential = np.zeros((10, 10))
potential[:, 0] = 1

potential = lib.laplace(potential, tolerance=1e-6)

fig, ax = plt.subplots(figsize=(10, 10), subplot_kw={"projection": "3d"})
X, Y = np.meshgrid(x, y)
surf = ax.plot_surface(X, Y, potential, cmap=cm.viridis)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('V')
plt.show()

'''
output is in the Q3.png File
'''
