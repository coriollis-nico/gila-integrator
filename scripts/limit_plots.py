#!/usr/bin/python

"""
Plots limit curves & finds point where dy/dx = 0
"""

# Modules
from os import makedirs
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('grayscale')
plt.rcParams['text.usetex'] = True

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

# Finds x value where dy/dx = 0
found_17 = False
found_22 = False
found_27 = False
for i in range(len(x[:-1])):
    if not found_17 and l17[i+1] == l17[i]:
        print("l17 transition point: a/a0 = {} in x[{}] and H/H0 = {}"
              .format(np.exp(x[i]), i, l17[i]))
        print("ln(a/a0) = {}".format(x[i]))
        x17 = x[i]
        found_17 = True
    if not found_22 and l22[i+1] == l22[i]:
        print("l22 transition point: a/a0 = {} in x[{}] and H/H0 = {}"
              .format(np.exp(x[i]), i, l22[i]))
        print("ln(a/a0) = {}".format(x[i]))
        x22 = x[i]
        found_22 = True
    if not found_27 and l27[i+1] == l27[i]:
        print("l27 transition point: a/a0 = {} in x[{}] and H/H0 = {}"
              .format(np.exp(x[i]), i, l27[i]))
        print("ln(a/a0) = {}".format(x[i]))
        x27 = x[i]
        found_27 = True

# Plotting
print("Plotting...")

plt.figure(layout="constrained")

plt.xscale("log")
plt.yscale("log")

plt.xlim(np.exp(x[-1]), np.exp(0))
plt.ylim(top=1.e33)

plt.plot(np.exp(x), l17, label=r"$ l = 1 \times 10^{-17} $",
         linestyle="dotted", color="k")
plt.plot(np.exp(x), l22, label=r"$ l = 1 \times 10^{-22} $",
         linestyle="dashed", color="k")
plt.plot(np.exp(x), l27, label=r"$ l = 1 \times 10^{-27} $",
         linestyle="dashdot", color="k")
plt.plot(np.exp(gr_integration[:, 0]), gr_integration[:, 1],
         linestyle="solid", label="RG", color="k")

plt.xlabel(r"$ \frac{a}{a_0} $")
plt.ylabel(r"$ \frac{H}{H_0} $")

plt.legend(loc="lower left", frameon=False)

plt.savefig(fig_dir+"/limits.pdf")

plt.close()

exit()
