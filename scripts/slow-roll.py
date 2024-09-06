#!/usr/bin/python

"""
Plots slow-roll data
"""

# Modules
from os import makedirs
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('grayscale')

# Data import
fig_dir = "plots/slowroll_plots"
data_dir = "data/sims/curve_sandwich"
makedirs(fig_dir, exist_ok=True)


sr_m03p01l17 = np.loadtxt(data_dir+"/slow-m03p01l1.0E-17.dat", comments='#')
sr_m03p01l22 = np.loadtxt(data_dir+"/slow-m03p01l1.0E-22.dat", comments='#')
sr_m03p01l27 = np.loadtxt(data_dir+"/slow-m03p01l1.0E-27.dat", comments='#')

sr_m08p02l17 = np.loadtxt(data_dir+"/slow-m08p02l1.0E-17.dat", comments='#')
sr_m08p02l22 = np.loadtxt(data_dir+"/slow-m08p02l1.0E-22.dat", comments='#')
sr_m08p02l27 = np.loadtxt(data_dir+"/slow-m08p02l1.0E-27.dat", comments='#')

sr_m10p10l17 = np.loadtxt(data_dir+"/slow-m10p10l1.0E-17.dat", comments='#')
sr_m10p10l22 = np.loadtxt(data_dir+"/slow-m10p10l1.0E-22.dat", comments='#')
sr_m10p10l27 = np.loadtxt(data_dir+"/slow-m10p10l1.0E-27.dat", comments='#')

l17_index = [i for i in range(len(sr_m03p01l17[:, 0]))
             if -60 <= sr_m03p01l17[i, 0] <= -50]
l22_index = [i for i in range(len(sr_m03p01l22[:, 0]))
             if -60 <= sr_m03p01l22[i, 0] <= -50]
l27_index = [i for i in range(len(sr_m03p01l27[:, 0]))
             if -60 <= sr_m03p01l27[i, 0] <= -50]

# Plotting
print("Plotting...")

fig, axs = plt.subplots(1, 3, sharey=True, sharex=True, layout="constrained",
                        figsize=[10., 3.5], dpi=400)

axs[0].set_title(r"$l = 1 \times 10^{-17}$")
axs[0].plot(sr_m03p01l17[l17_index, 0], sr_m03p01l17[l17_index, 2],
            label=r"$ m=3, p=1 $", color="k", linestyle="solid")
axs[0].plot(sr_m08p02l17[l17_index, 0], sr_m08p02l17[l17_index, 2],
            label=r"$ m=8, p=2 $", color="k", linestyle="dashed")
axs[0].plot(sr_m10p10l17[l17_index, 0], sr_m10p10l17[l17_index, 2],
            label=r"$ m=10, p=10 $", color="k", linestyle="dotted")

axs[1].set_title(r"$l = 1 \times 10^{-22}$")
axs[1].plot(sr_m03p01l22[l22_index, 0], sr_m03p01l22[l22_index, 2],
            label=r"$ m=3, p=1 $", color="k", linestyle="solid")
axs[1].plot(sr_m08p02l22[l22_index, 0], sr_m08p02l22[l22_index, 2],
            label=r"$ m=8, p=2 $", color="k", linestyle="dashed")
axs[1].plot(sr_m10p10l22[l22_index, 0], sr_m10p10l22[l22_index, 2],
            label=r"$ m=10, p=10 $", color="k", linestyle="dotted")

axs[2].set_title(r"$l = 1 \times 10^{-27}$")
axs[2].plot(sr_m03p01l27[l27_index, 0], sr_m03p01l27[l27_index, 2],
            label=r"$ m=3, p=1 $", color="k", linestyle="solid")
axs[2].plot(sr_m08p02l27[l27_index, 0], sr_m08p02l27[l27_index, 2],
            label=r"$ m=8, p=2 $", color="k", linestyle="dashed")
axs[2].plot(sr_m10p10l27[l27_index, 0], sr_m10p10l27[l27_index, 2],
            label=r"$ m=10, p=10 $", color="k", linestyle="dotted")

axs[1].set_xlabel(r"$N$")
axs[0].set_ylabel(r"$ \epsilon_1 $")
for ax in axs:
    ax.set_xlim(-60, -50)
    ax.set_ylim(0)
    ax.legend(loc="best", fontsize=8)
    ax.fill_between([-60, -50], 0.0097, -0.01, alpha=0.08)

plt.savefig(fig_dir+"/e1.png")

plt.close()

fig, axs = plt.subplots(1, 3, sharey=True, sharex=True, layout="constrained",
                        figsize=[10., 3.5], dpi=400)

axs[0].set_title(r"$l = 1 \times 10^{-17}$")
axs[0].plot(sr_m03p01l17[l17_index, 0], sr_m03p01l17[l17_index, 3],
            label=r"$ m=3, p=1 $", color="k", linestyle="solid")
axs[0].plot(sr_m08p02l17[l17_index, 0], sr_m08p02l17[l17_index, 3],
            label=r"$ m=8, p=2 $", color="k", linestyle="dashed")
axs[0].plot(sr_m10p10l17[l17_index, 0], sr_m10p10l17[l17_index, 3],
            label=r"$ m=10, p=10 $", color="k", linestyle="dotted")

axs[1].set_title(r"$l = 1 \times 10^{-22}$")
axs[1].plot(sr_m03p01l22[l22_index, 0], sr_m03p01l22[l22_index, 3],
            label=r"$ m=3, p=1 $", color="k", linestyle="solid")
axs[1].plot(sr_m08p02l22[l22_index, 0], sr_m08p02l22[l22_index, 3],
            label=r"$ m=8, p=2 $", color="k", linestyle="dashed")
axs[1].plot(sr_m10p10l22[l22_index, 0], sr_m10p10l22[l22_index, 3],
            label=r"$ m=10, p=10 $", color="k", linestyle="dotted")

axs[2].set_title(r"$l = 1 \times 10^{-27}$")
axs[2].plot(sr_m03p01l27[l27_index, 0], sr_m03p01l27[l27_index, 3],
            label=r"$ m=3, p=1 $", color="k", linestyle="solid")
axs[2].plot(sr_m08p02l27[l27_index, 0], sr_m08p02l27[l27_index, 3],
            label=r"$ m=8, p=2 $", color="k", linestyle="dashed")
axs[2].plot(sr_m10p10l27[l27_index, 0], sr_m10p10l27[l27_index, 3],
            label=r"$ m=10, p=10 $", color="k", linestyle="dotted")

axs[1].set_xlabel(r"$N$")
axs[0].set_ylabel(r"$ \epsilon_2 $")
for ax in axs:
    ax.set_xlim(-60, -50)
    ax.set_ylim(0.016, 0.0215)
    ax.legend(loc="best", fontsize=8)
    ax.fill_between([-60, -50], 0.032+0.009, 0.032-0.008, alpha=0.05)

plt.savefig(fig_dir+"/e2.png")

plt.close()

fig, axs = plt.subplots(1, 3, sharey=True, sharex=True, layout="constrained",
                        figsize=[10., 3.5], dpi=400)

axs[0].set_title(r"$l = 1 \times 10^{-17}$")
axs[0].plot(sr_m03p01l17[l17_index, 0], sr_m03p01l17[l17_index, 4],
            label=r"$ m=3, p=1 $", color="k", linestyle="solid")
axs[0].plot(sr_m08p02l17[l17_index, 0], sr_m08p02l17[l17_index, 4],
            label=r"$ m=8, p=2 $", color="k", linestyle="dashed")
axs[0].plot(sr_m10p10l17[l17_index, 0], sr_m10p10l17[l17_index, 4],
            label=r"$ m=10, p=10 $", color="k", linestyle="dotted")

axs[1].set_title(r"$l = 1 \times 10^{-22}$")
axs[1].plot(sr_m03p01l22[l22_index, 0], sr_m03p01l22[l22_index, 4],
            label=r"$ m=3, p=1 $", color="k", linestyle="solid")
axs[1].plot(sr_m08p02l22[l22_index, 0], sr_m08p02l22[l22_index, 4],
            label=r"$ m=8, p=2 $", color="k", linestyle="dashed")
axs[1].plot(sr_m10p10l22[l22_index, 0], sr_m10p10l22[l22_index, 4],
            label=r"$ m=10, p=10 $", color="k", linestyle="dotted")

axs[2].set_title(r"$l = 1 \times 10^{-27}$")
axs[2].plot(sr_m03p01l27[l27_index, 0], sr_m03p01l27[l27_index, 4],
            label=r"$ m=3, p=1 $", color="k", linestyle="solid")
axs[2].plot(sr_m08p02l27[l27_index, 0], sr_m08p02l27[l27_index, 4],
            label=r"$ m=8, p=2 $", color="k", linestyle="dashed")
axs[2].plot(sr_m10p10l27[l27_index, 0], sr_m10p10l27[l27_index, 4],
            label=r"$ m=10, p=10 $", color="k", linestyle="dotted")

axs[1].set_xlabel(r"$N$")
axs[0].set_ylabel(r"$ \epsilon_2 $")
for ax in axs:
    ax.set_xlim(-60, -50)
    ax.set_ylim(0.0165, 0.0215)
    ax.legend(loc="best", fontsize=8)
    ax.fill_between([-60, -50], 0.19+0.55, 0.19-0.53, alpha=0.08)

plt.savefig(fig_dir+"/e3.png")

plt.close()
