#!/usr/bin/python

"""
Plots solutions
"""

# Modules
from os import makedirs
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('grayscale')

# Data import
fig_dir = "plots/solution_plots"
data_dir = "data/sims/curve_sandwich"
makedirs(fig_dir, exist_ok=True)

print("Reading sandwich_curves")
x = np.loadtxt(data_dir+"/x_grid.dat", comments='#')

m03p01l17 = np.loadtxt(data_dir+"/m03p01l1.0E-17.dat", comments='#')
m03p01l22 = np.loadtxt(data_dir+"/m03p01l1.0E-22.dat", comments='#')
m03p01l27 = np.loadtxt(data_dir+"/m03p01l1.0E-27.dat", comments='#')

m08p02l17 = np.loadtxt(data_dir+"/m08p02l1.0E-17.dat", comments='#')
m08p02l22 = np.loadtxt(data_dir+"/m08p02l1.0E-22.dat", comments='#')
m08p02l27 = np.loadtxt(data_dir+"/m08p02l1.0E-27.dat", comments='#')

m10p10l17 = np.loadtxt(data_dir+"/m10p10l1.0E-17.dat", comments='#')
m10p10l22 = np.loadtxt(data_dir+"/m10p10l1.0E-22.dat", comments='#')
m10p10l27 = np.loadtxt(data_dir+"/m10p10l1.0E-27.dat", comments='#')

# Plotting
print("Plotting...")

fig, axs = plt.subplots(1, 3, sharey=True, sharex=True, layout="constrained",
                        figsize=[10., 3.5], dpi=400)


axs[0].set_title(r"$l = 1 \times 10^{-17}$")
axs[0].plot(np.exp(x), m03p01l17, label=r"$m=3, p=1$")
axs[0].plot(np.exp(x), m08p02l17, label=r"$m=8, p=2$")
axs[0].plot(np.exp(x), m10p10l17, label=r"$m=10, p=10$")

axs[1].set_title(r"$l = 1 \times 10^{-22}$")
axs[1].plot(np.exp(x), m03p01l22, label=r"$m=3, p=1$")
axs[1].plot(np.exp(x), m08p02l22, label=r"$m=8, p=2$")
axs[1].plot(np.exp(x), m10p10l22, label=r"$m=10, p=10$")

axs[2].set_title(r"$l = 1 \times 10^{-27}$")
axs[2].plot(np.exp(x), m03p01l27, label=r"$m=3, p=1$")
axs[2].plot(np.exp(x), m08p02l27, label=r"$m=8, p=2$")
axs[2].plot(np.exp(x), m10p10l27, label=r"$m=10, p=10$")

plt.yscale("log")
axs[1].set_xlabel(r"$\frac{a}{a_0}$")
axs[0].set_ylabel(r"$ \frac{H}{H_0} $")
for ax in axs:
    ax.set_xscale("log")
    ax.set_xlim(np.exp(x[-1]), np.exp(x[0]))
    ax.legend(loc="lower left", frameon=False)

plt.savefig(fig_dir+"/solutions.png")

plt.close()

exit()
