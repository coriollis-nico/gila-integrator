#!/usr/bin/python

"""
Plots extended slow-roll data
"""

# Modules
from os import makedirs
import numpy as np
import matplotlib.pyplot as plt

plt.style.use("grayscale")
plt.rcParams["text.usetex"] = True
plt.rcParams["figure.figsize"] = [6.4, 2.4]

# Data import
fig_dir = "plots/extended_slowroll_plots"
data_dir = "data/sims/curve_sandwich"
makedirs(fig_dir, exist_ok=True)

lim_low = -60
lim_up = 0


# Planck data
e1_max = 0.0097
e2 = 0.032
e2_min = e2 - 0.008
e2_max = e2 + 0.009
e3 = 0.19
e3_min = e3 - 0.53
e3_max = e3 + 0.55


sr_m03p01l17 = np.loadtxt(data_dir + "/slow-m03p01l1.0E-17.dat", comments="#")
sr_m03p01l22 = np.loadtxt(data_dir + "/slow-m03p01l1.0E-22.dat", comments="#")
sr_m03p01l27 = np.loadtxt(data_dir + "/slow-m03p01l1.0E-27.dat", comments="#")

sr_m08p02l17 = np.loadtxt(data_dir + "/slow-m08p02l1.0E-17.dat", comments="#")
sr_m08p02l22 = np.loadtxt(data_dir + "/slow-m08p02l1.0E-22.dat", comments="#")
sr_m08p02l27 = np.loadtxt(data_dir + "/slow-m08p02l1.0E-27.dat", comments="#")

sr_m10p10l17 = np.loadtxt(data_dir + "/slow-m10p10l1.0E-17.dat", comments="#")
sr_m10p10l22 = np.loadtxt(data_dir + "/slow-m10p10l1.0E-22.dat", comments="#")
sr_m10p10l27 = np.loadtxt(data_dir + "/slow-m10p10l1.0E-27.dat", comments="#")

sr_m99p99l17 = np.loadtxt(data_dir + "/slow-m99p99l1.0E-17.dat", comments="#")
sr_m99p99l22 = np.loadtxt(data_dir + "/slow-m99p99l1.0E-22.dat", comments="#")
sr_m99p99l27 = np.loadtxt(data_dir + "/slow-m99p99l1.0E-27.dat", comments="#")

l17_index = [
    i for i in range(len(sr_m03p01l17[:, 0])) if lim_low <= sr_m03p01l17[i, 0] <= lim_up
]
l22_index = [
    i for i in range(len(sr_m03p01l22[:, 0])) if lim_low <= sr_m03p01l22[i, 0] <= lim_up
]
l27_index = [
    i for i in range(len(sr_m03p01l27[:, 0])) if lim_low <= sr_m03p01l27[i, 0] <= lim_up
]

# Plotting
print("Plotting...")

fig, axs = plt.subplots(1, 3, sharey=True, sharex=True, layout="constrained")

axs[0].set_title(r"$l = 1 \times 10^{-17}$")
axs[0].plot(
    sr_m03p01l17[l17_index, 0],
    sr_m03p01l17[l17_index, 2],
    label=r"$ m=3, p=1 $",
    color="k",
    linestyle="solid",
)
axs[0].plot(
    sr_m08p02l17[l17_index, 0],
    sr_m08p02l17[l17_index, 2],
    label=r"$ m=8, p=2 $",
    color="k",
    linestyle="dashed",
)
axs[0].plot(
    sr_m10p10l17[l17_index, 0],
    sr_m10p10l17[l17_index, 2],
    label=r"$ m=10, p=10 $",
    color="k",
    linestyle="dotted",
)
axs[0].plot(
    sr_m99p99l17[l17_index, 0],
    sr_m99p99l17[l17_index, 2],
    label=r"$ m=99, p=99 $",
    color="k",
    linestyle="dashdot",
)

axs[1].set_title(r"$l = 1 \times 10^{-22}$")
axs[1].plot(
    sr_m03p01l22[l22_index, 0],
    sr_m03p01l22[l22_index, 2],
    label=r"$ m=3, p=1 $",
    color="k",
    linestyle="solid",
)
axs[1].plot(
    sr_m08p02l22[l22_index, 0],
    sr_m08p02l22[l22_index, 2],
    label=r"$ m=8, p=2 $",
    color="k",
    linestyle="dashed",
)
axs[1].plot(
    sr_m10p10l22[l22_index, 0],
    sr_m10p10l22[l22_index, 2],
    label=r"$ m=10, p=10 $",
    color="k",
    linestyle="dotted",
)
axs[1].plot(
    sr_m99p99l22[l22_index, 0],
    sr_m99p99l22[l22_index, 2],
    label=r"$ m=99, p=99 $",
    color="k",
    linestyle="dashdot",
)

axs[2].set_title(r"$l = 1 \times 10^{-27}$")
axs[2].plot(
    sr_m03p01l27[l27_index, 0],
    sr_m03p01l27[l27_index, 2],
    label=r"$ m=3, p=1 $",
    color="k",
    linestyle="solid",
)
axs[2].plot(
    sr_m08p02l27[l27_index, 0],
    sr_m08p02l27[l27_index, 2],
    label=r"$ m=8, p=2 $",
    color="k",
    linestyle="dashed",
)
axs[2].plot(
    sr_m10p10l27[l27_index, 0],
    sr_m10p10l27[l27_index, 2],
    label=r"$ m=10, p=10 $",
    color="k",
    linestyle="dotted",
)
axs[2].plot(
    sr_m99p99l27[l27_index, 0],
    sr_m99p99l27[l27_index, 2],
    label=r"$ m=99, p=99 $",
    color="k",
    linestyle="dashdot",
)

# axs[1].set_xlabel(r"$N$")
axs[0].set_ylabel(r"$ \epsilon_1 $")
axs[0].legend(loc="upper left", fontsize=8)
for ax in axs:
    ax.set_xlim(-0.005, lim_up)
    ax.set_ylim(bottom=0)
    ax.fill_between([lim_low, lim_up], e1_max, -1, alpha=0.14)

plt.savefig(fig_dir + "/e1.pdf")

plt.close()

fig, axs = plt.subplots(1, 3, sharey=True, sharex=True, layout="constrained")

# axs[0].set_title(r"$l = 1 \times 10^{-17}$")
axs[0].plot(
    sr_m03p01l17[l17_index, 0],
    sr_m03p01l17[l17_index, 3],
    label=r"$ m=3, p=1 $",
    color="k",
    linestyle="solid",
)
axs[0].plot(
    sr_m08p02l17[l17_index, 0],
    sr_m08p02l17[l17_index, 3],
    label=r"$ m=8, p=2 $",
    color="k",
    linestyle="dashed",
)
axs[0].plot(
    sr_m10p10l17[l17_index, 0],
    sr_m10p10l17[l17_index, 3],
    label=r"$ m=10, p=10 $",
    color="k",
    linestyle="dotted",
)
axs[0].plot(
    sr_m99p99l17[l17_index, 0],
    sr_m99p99l17[l17_index, 3],
    label=r"$ m=99, p=99 $",
    color="k",
    linestyle="dashdot",
)

# axs[1].set_title(r"$l = 1 \times 10^{-22}$")
axs[1].plot(
    sr_m03p01l22[l22_index, 0],
    sr_m03p01l22[l22_index, 3],
    label=r"$ m=3, p=1 $",
    color="k",
    linestyle="solid",
)
axs[1].plot(
    sr_m08p02l22[l22_index, 0],
    sr_m08p02l22[l22_index, 3],
    label=r"$ m=8, p=2 $",
    color="k",
    linestyle="dashed",
)
axs[1].plot(
    sr_m10p10l22[l22_index, 0],
    sr_m10p10l22[l22_index, 3],
    label=r"$ m=10, p=10 $",
    color="k",
    linestyle="dotted",
)
axs[1].plot(
    sr_m99p99l22[l22_index, 0],
    sr_m99p99l22[l22_index, 3],
    label=r"$ m=99, p=99 $",
    color="k",
    linestyle="dashdot",
)

# axs[2].set_title(r"$l = 1 \times 10^{-27}$")
axs[2].plot(
    sr_m03p01l27[l27_index, 0],
    sr_m03p01l27[l27_index, 3],
    label=r"$ m=3, p=1 $",
    color="k",
    linestyle="solid",
)
axs[2].plot(
    sr_m08p02l27[l27_index, 0],
    sr_m08p02l27[l27_index, 3],
    label=r"$ m=8, p=2 $",
    color="k",
    linestyle="dashed",
)
axs[2].plot(
    sr_m10p10l27[l27_index, 0],
    sr_m10p10l27[l27_index, 3],
    label=r"$ m=10, p=10 $",
    color="k",
    linestyle="dotted",
)
axs[2].plot(
    sr_m99p99l27[l27_index, 0],
    sr_m99p99l27[l27_index, 3],
    label=r"$ m=99, p=99 $",
    color="k",
    linestyle="dashdot",
)

# axs[1].set_xlabel(r"$N$")
axs[0].set_ylabel(r"$ \epsilon_2 $")
for ax in axs:
    ax.set_xlim(lim_low, lim_up)
    # ax.legend(loc="best", fontsize=8)
    ax.fill_between([lim_low, lim_up], e2_min, e2_max, alpha=0.1)

plt.savefig(fig_dir + "/e2.pdf")

plt.close()

fig, axs = plt.subplots(1, 3, sharey=True, sharex=True, layout="constrained")

# axs[0].set_title(r"$l = 1 \times 10^{-17}$")
axs[0].plot(
    sr_m03p01l17[l17_index, 0],
    sr_m03p01l17[l17_index, 4],
    label=r"$ m=3, p=1 $",
    color="k",
    linestyle="solid",
)
axs[0].plot(
    sr_m08p02l17[l17_index, 0],
    sr_m08p02l17[l17_index, 4],
    label=r"$ m=8, p=2 $",
    color="k",
    linestyle="dashed",
)
axs[0].plot(
    sr_m10p10l17[l17_index, 0],
    sr_m10p10l17[l17_index, 4],
    label=r"$ m=10, p=10 $",
    color="k",
    linestyle="dotted",
)
axs[0].plot(
    sr_m99p99l17[l17_index, 0],
    sr_m99p99l17[l17_index, 4],
    label=r"$ m=99, p=99 $",
    color="k",
    linestyle="dashdot",
)

# axs[1].set_title(r"$l = 1 \times 10^{-22}$")
axs[1].plot(
    sr_m03p01l22[l22_index, 0],
    sr_m03p01l22[l22_index, 4],
    label=r"$ m=3, p=1 $",
    color="k",
    linestyle="solid",
)
axs[1].plot(
    sr_m08p02l22[l22_index, 0],
    sr_m08p02l22[l22_index, 4],
    label=r"$ m=8, p=2 $",
    color="k",
    linestyle="dashed",
)
axs[1].plot(
    sr_m10p10l22[l22_index, 0],
    sr_m10p10l22[l22_index, 4],
    label=r"$ m=10, p=10 $",
    color="k",
    linestyle="dotted",
)
axs[1].plot(
    sr_m99p99l22[l22_index, 0],
    sr_m99p99l22[l22_index, 4],
    label=r"$ m=99, p=99 $",
    color="k",
    linestyle="dashdot",
)

# axs[2].set_title(r"$l = 1 \times 10^{-27}$")
axs[2].plot(
    sr_m03p01l27[l27_index, 0],
    sr_m03p01l27[l27_index, 4],
    label=r"$ m=3, p=1 $",
    color="k",
    linestyle="solid",
)
axs[2].plot(
    sr_m08p02l27[l27_index, 0],
    sr_m08p02l27[l27_index, 4],
    label=r"$ m=8, p=2 $",
    color="k",
    linestyle="dashed",
)
axs[2].plot(
    sr_m10p10l27[l27_index, 0],
    sr_m10p10l27[l27_index, 4],
    label=r"$ m=10, p=10 $",
    color="k",
    linestyle="dotted",
)
axs[2].plot(
    sr_m99p99l27[l27_index, 0],
    sr_m99p99l27[l27_index, 4],
    label=r"$ m=99, p=99 $",
    color="k",
    linestyle="dashdot",
)

axs[1].set_xlabel(r"$N$")
axs[0].set_ylabel(r"$ \epsilon_3 $")
for ax in axs:
    ax.set_xlim(lim_low, lim_up)
    # ax.legend(loc="best", fontsize=8)
    ax.fill_between([lim_low, lim_up], e3_min, e3_max, alpha=0.14)

plt.savefig(fig_dir + "/e3.pdf")

plt.close()
