import numpy as np
import matplotlib.pyplot as plt
import lib


nx = 20+1  # number of points in the length dimension
nt = 5000+1  # number of points in the time dimension
Lx = 2    # domain length in the length dimension
dt = 0.008  # grid spacing in the time dimension
dx = Lx/(nx-1)  # grid spacing in the length dimension
Lt = (nt-1)*dt    # domain length in the time dimension

x = np.linspace(0, Lx, nx)
t = np.linspace(0, Lt, nt)
alpha = dt/(2*dx**2)
boundary_values = [0, 0]
boundary_time = 500

u = np.zeros((nt, nx))  # array to store the temperature over time

u[0, :] = 20*np.abs(np.sin(np.pi*x))

heat = lib.heat_explicit(u.copy(), nt, alpha, boundary_values, boundary_time)

fig, ax = plt.subplots(figsize=(12, 8))
times = [0, 10, 20, 50, 100, 200, 500]
fig, axs = plt.subplots(nrows=4, ncols=3, figsize=(15, 12))
plt.subplots_adjust(hspace=0.5)
fig.suptitle("Residue versus iteration steps", fontsize=18, y=0.95)

# axs[1][1].plot(x, heat[0, :])
# axs[1][1].set_title("t=0")
# axs[1][1].set_ylabel('T')
# axs[1][1].set_xlabel('X')
#
# axs[1][2].plot(x, heat[10, :])
# axs[1][2].set_title("t=0")
# axs[1][2].set_ylabel('T')
# axs[1][2].set_xlabel('X')
#
# axs[2][1].plot(x, heat[20, :])
# axs[2][1].set_title("t=0")
# axs[2][1].set_ylabel('T')
# axs[2][1].set_xlabel('X')
#
# axs[4].plot(x, heat[50, :])
# axs[4].set_title("t=0")
# axs[4].set_ylabel('T')
# axs[4].set_xlabel('X')

plt.show()
