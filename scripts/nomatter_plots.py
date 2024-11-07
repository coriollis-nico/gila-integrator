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
fig_dir = "plots/nomatter_plots"
data_dir = "data/sims/no_matter"
makedirs(fig_dir, exist_ok=True)

print("Reading data")
x = np.loadtxt(data_dir+"/x_grid.dat", comments='#')

m03p01 = np.loadtxt(data_dir+"/m03p01.dat", comments='#')
m08p02 = np.loadtxt(data_dir+"/m08p02.dat", comments='#')
m10p10 = np.loadtxt(data_dir+"/m10p10.dat", comments='#')

# Plotting
print("Plotting...")

plt.figure(layout="constrained", dpi=400)

plt.plot(np.exp(x), m03p01, label=r"$m=3, p=1$", color="k",
         linestyle="solid")
plt.plot(np.exp(x), m08p02, label=r"$m=8, p=2$", color="k",
         linestyle="dashed")
plt.plot(np.exp(x), m08p02, label=r"$m=10, p=10$", color="k",
         linestyle="dotted")

plt.yscale("log")
plt.xscale("log")

plt.xlabel(r"$\frac{a}{a_i}$")
plt.ylabel(r"$ \frac{H}{H_i} $")

plt.xlim(np.exp(x[-1]), np.exp(x[0]))
plt.legend(loc="lower left", frameon=False)

plt.savefig(fig_dir+"/solutions.png")

plt.close()

exit()
