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
fig_dir = "plots/limit_plots"
data_dir = "data/sims/curve_sandwich"
makedirs(fig_dir, exist_ok=True)


print("Reading gr.dat")
gr_integration = np.loadtxt("data/sims/gr/gr.dat", comments='#')

print("Reading sandwich_curves")
x = np.loadtxt(data_dir+"/x_grid.dat", comments='#')

l17 = np.loadtxt(data_dir+"/lim_l1.0E-17.dat", comments='#')
l22 = np.loadtxt(data_dir+"/lim_l1.0E-22.dat", comments='#')
l27 = np.loadtxt(data_dir+"/lim_l1.0E-27.dat", comments='#')

# Plotting
print("Plotting...")

plt.figure(dpi=400, layout="constrained")

plt.xscale("log")
plt.yscale("log")

plt.xlim(np.exp(x[-1]), np.exp(0))

plt.plot(np.exp(x), l17, label=r"$ l = 1 \times 10^{-17} $")
plt.plot(np.exp(x), l22, label=r"$ l = 1 \times 10^{-22} $")
plt.plot(np.exp(x), l27, label=r"$ l = 1 \times 10^{-27} $")
plt.plot(np.exp(gr_integration[:, 0]), gr_integration[:, 1], "--", label="RG")

plt.xlabel(r"$ \frac{a}{a_0} $")
plt.ylabel(r"$ \frac{H}{H_0} $")

plt.legend(loc="lower left", frameon=False)

plt.savefig(fig_dir+"/limits.png")

plt.close()

exit()
