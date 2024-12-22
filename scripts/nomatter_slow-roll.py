#!/usr/bin/python

"""
Plots no_matter slow-roll data
"""

# Modules
from os import makedirs
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('grayscale')
plt.rcParams['text.usetex'] = True

# Data import
fig_dir = "plots/nomatter_plots"
data_dir = "data/sims/no_matter"
makedirs(fig_dir, exist_ok=True)


sr_m03p01 = np.loadtxt(data_dir+"/slow-m03p01.dat", comments='#')
sr_m08p02 = np.loadtxt(data_dir+"/slow-m08p02.dat", comments='#')
sr_m10p10 = np.loadtxt(data_dir+"/slow-m10p10.dat", comments='#')

# Plotting
print("Plotting...")

# ϵ_1
plt.figure(layout="constrained", dpi=250)

plt.plot(sr_m03p01[:, 0], sr_m03p01[:, 2],
         label=r"$ m=3, p=1 $", color="k", linestyle="solid")
plt.plot(sr_m08p02[:, 0], sr_m08p02[:, 2],
         label=r"$ m=8, p=2 $", color="k", linestyle="dashed")
plt.plot(sr_m10p10[:, 0], sr_m10p10[:, 2],
         label=r"$ m=3, p=1 $", color="k", linestyle="dotted")

plt.xlabel(r"$N$")
plt.ylabel(r"$ \epsilon_1 $")
plt.legend(loc="best", fontsize=8)
plt.fill_between([-60, 0], 0.0097, -0.01, alpha=0.16)
plt.ylim(-0.01)

plt.savefig(fig_dir+"/e1.png")

plt.close()

# ϵ_2

plt.figure(layout="constrained", dpi=250)

plt.plot(sr_m03p01[:, 0], sr_m03p01[:, 3],
         label=r"$ m=3, p=1 $", color="k", linestyle="solid")
plt.plot(sr_m08p02[:, 0], sr_m08p02[:, 3],
         label=r"$ m=8, p=2 $", color="k", linestyle="dashed")
plt.plot(sr_m10p10[:, 0], sr_m10p10[:, 3],
         label=r"$ m=3, p=1 $", color="k", linestyle="dotted")

plt.xlabel(r"$N$")
plt.ylabel(r"$ \epsilon_2 $")
plt.legend(loc="best", fontsize=8)
plt.fill_between([-60, 0], 0.032+0.009, 0.032-0.008, alpha=0.1)

plt.savefig(fig_dir+"/e2.png")

plt.close()

# ϵ_3
plt.figure(layout="constrained", dpi=250)

plt.plot(sr_m03p01[:, 0], sr_m03p01[:, 4],
         label=r"$ m=3, p=1 $", color="k", linestyle="solid")
plt.plot(sr_m08p02[:, 0], sr_m08p02[:, 4],
         label=r"$ m=8, p=2 $", color="k", linestyle="dashed")
plt.plot(sr_m10p10[:, 0], sr_m10p10[:, 4],
         label=r"$ m=3, p=1 $", color="k", linestyle="dotted")

plt.xlabel(r"$N$")
plt.ylabel(r"$ \epsilon_3 $")
plt.legend(loc="best", fontsize=8)
plt.fill_between([-60, 0], 0.19+0.55, 0.19-0.53, alpha=0.16)

plt.savefig(fig_dir+"/e3.png")

plt.close()
