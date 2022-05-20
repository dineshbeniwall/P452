import numpy as np
import matplotlib.pyplot as plt
import lib


nx = 20+1
nt = 5000+1
Lx = 2
dt = 0.008
dx = Lx/(nx-1)
Lt = (nt-1)*dt

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
fig, axs = plt.subplots(nrows=4, ncols=2, figsize=(15, 12))
plt.subplots_adjust(hspace=1)
fig.suptitle("temperature at diff. time", fontsize=18, y=0.95)

axs[0][0].plot(x, heat[0, :])
axs[0][0].set_title("t=0")
axs[0][0].set_ylabel('T')
axs[0][0].set_xlabel('X')

axs[0][1].plot(x, heat[10, :])
axs[0][1].set_title("t=10")
axs[0][1].set_ylabel('T')
axs[0][1].set_xlabel('X')
#
axs[1][0].plot(x, heat[20, :])
axs[1][0].set_title("t=20")
axs[1][0].set_ylabel('T')
axs[1][0].set_xlabel('X')
#
axs[1][1].plot(x, heat[50, :])
axs[1][1].set_title("t=50")
axs[1][1].set_ylabel('T')
axs[1][1].set_xlabel('X')

axs[2][0].plot(x, heat[100, :])
axs[2][0].set_title("t=100")
axs[2][0].set_ylabel('T')
axs[2][0].set_xlabel('X')
#
axs[2][1].plot(x, heat[200, :])
axs[2][1].set_title("t=200")
axs[2][1].set_ylabel('T')
axs[2][1].set_xlabel('X')

axs[3][0].plot(x, heat[500, :])
axs[3][0].set_title("t=500")
axs[3][0].set_ylabel('T')
axs[3][0].set_xlabel('X')

plt.show()
