#!/usr/bin/python

"""
Plots slow-roll data
"""

# Modules
from os import makedirs
import numpy as np
import matplotlib.pyplot as plt

plt.style.use("grayscale")
plt.rcParams["text.usetex"] = True
plt.rcParams["figure.figsize"] = [6.4, 3]

# Data import
fig_dir = "plots/slowroll_extra_plots"
data_dir = "data/sims/extra_cases"
makedirs(fig_dir, exist_ok=True)

# Planck data
e1_max = 0.0097
e2 = 0.032
e2_min = e2 - 0.008
e2_max = e2 + 0.009
e3 = 0.19
e3_min = e3 - 0.53
e3_max = e3 + 0.55


sr_m87p05l17 = np.loadtxt(data_dir + "/slow-m87p05l1.0E-17.dat", comments="#")
sr_m87p05l22 = np.loadtxt(data_dir + "/slow-m87p05l1.0E-22.dat", comments="#")
sr_m87p05l27 = np.loadtxt(data_dir + "/slow-m87p05l1.0E-27.dat", comments="#")

sr_m88p06l17 = np.loadtxt(data_dir + "/slow-m88p06l1.0E-17.dat", comments="#")
sr_m88p06l22 = np.loadtxt(data_dir + "/slow-m88p06l1.0E-22.dat", comments="#")
sr_m88p06l27 = np.loadtxt(data_dir + "/slow-m88p06l1.0E-27.dat", comments="#")

sr_m89p07l17 = np.loadtxt(data_dir + "/slow-m89p07l1.0E-17.dat", comments="#")
sr_m89p07l22 = np.loadtxt(data_dir + "/slow-m89p07l1.0E-22.dat", comments="#")
sr_m89p07l27 = np.loadtxt(data_dir + "/slow-m89p07l1.0E-27.dat", comments="#")

l17_index = [
    i for i in range(len(sr_m87p05l17[:, 0])) if -60 <= sr_m87p05l17[i, 0] <= -50
]
l22_index = [
    i for i in range(len(sr_m87p05l22[:, 0])) if -60 <= sr_m87p05l22[i, 0] <= -50
]
l27_index = [
    i for i in range(len(sr_m87p05l27[:, 0])) if -60 <= sr_m87p05l27[i, 0] <= -50
]

# Plotting
print("Plotting...")

fig, axs = plt.subplots(1, 3, sharey=True, sharex=True, layout="constrained")

axs[0].set_title(r"$l = 1 \times 10^{-17}$")
axs[0].plot(
    sr_m87p05l17[l17_index, 0],
    sr_m87p05l17[l17_index, 2],
    label=r"$ m=87, p=5   $",
    color="k",
    linestyle="solid",
)
axs[0].plot(
    sr_m88p06l17[l17_index, 0],
    sr_m88p06l17[l17_index, 2],
    label=r"$ m=88, p=6   $",
    color="k",
    linestyle="dashed",
)
axs[0].plot(
    sr_m89p07l17[l17_index, 0],
    sr_m89p07l17[l17_index, 2],
    label=r"$ m=89, p=7 $",
    color="k",
    linestyle="dotted",
)

axs[1].set_title(r"$l = 1 \times 10^{-22}$")
axs[1].plot(
    sr_m87p05l22[l22_index, 0],
    sr_m87p05l22[l22_index, 2],
    label=r"$ m=87, p=5 $",
    color="k",
    linestyle="solid",
)
axs[1].plot(
    sr_m88p06l22[l22_index, 0],
    sr_m88p06l22[l22_index, 2],
    label=r"$ m=88, p=6 $",
    color="k",
    linestyle="dashed",
)
axs[1].plot(
    sr_m89p07l22[l22_index, 0],
    sr_m89p07l22[l22_index, 2],
    label=r"$ m=89, p=7 $",
    color="k",
    linestyle="dotted",
)

axs[2].set_title(r"$l = 1 \times 10^{-27}$")
axs[2].plot(
    sr_m87p05l27[l27_index, 0],
    sr_m87p05l27[l27_index, 2],
    label=r"$ m=87, p=5 $",
    color="k",
    linestyle="solid",
)
axs[2].plot(
    sr_m88p06l27[l27_index, 0],
    sr_m88p06l27[l27_index, 2],
    label=r"$ m=88, p=6 $",
    color="k",
    linestyle="dashed",
)
axs[2].plot(
    sr_m89p07l27[l27_index, 0],
    sr_m89p07l27[l27_index, 2],
    label=r"$ m=89, p=7 $",
    color="k",
    linestyle="dotted",
)

axs[1].set_xlabel(r"$N$")
axs[0].legend(loc="upper left", fontsize=7)
axs[0].set_ylabel(r"$ \epsilon_1 $")
for ax in axs:
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.set_xlim(-60, -50)
    ax.set_ylim(bottom=0)
    ax.fill_between([-60, -50], e1_max, 0, alpha=0.14)

plt.savefig(fig_dir + "/e1.pdf")

plt.close()

fig, axs = plt.subplots(1, 3, sharey=True, sharex=True, layout="constrained")

axs[0].set_title(r"$l = 1 \times 10^{-17}$")
axs[0].plot(
    sr_m87p05l17[l17_index, 0],
    sr_m87p05l17[l17_index, 3],
    label=r"$ m=87, p=5 $",
    color="k",
    linestyle="solid",
)
axs[0].plot(
    sr_m88p06l17[l17_index, 0],
    sr_m88p06l17[l17_index, 3],
    label=r"$ m=88, p=6 $",
    color="k",
    linestyle="dashed",
)
axs[0].plot(
    sr_m89p07l17[l17_index, 0],
    sr_m89p07l17[l17_index, 3],
    label=r"$ m=89, p=7 $",
    color="k",
    linestyle="dotted",
)

axs[1].set_title(r"$l = 1 \times 10^{-22}$")
axs[1].plot(
    sr_m87p05l22[l22_index, 0],
    sr_m87p05l22[l22_index, 3],
    label=r"$ m=87, p=5 $",
    color="k",
    linestyle="solid",
)
axs[1].plot(
    sr_m88p06l22[l22_index, 0],
    sr_m88p06l22[l22_index, 3],
    label=r"$ m=88, p=6 $",
    color="k",
    linestyle="dashed",
)
axs[1].plot(
    sr_m89p07l22[l22_index, 0],
    sr_m89p07l22[l22_index, 3],
    label=r"$ m=89, p=7 $",
    color="k",
    linestyle="dotted",
)

axs[2].set_title(r"$l = 1 \times 10^{-27}$")
axs[2].plot(
    sr_m87p05l27[l27_index, 0],
    sr_m87p05l27[l27_index, 3],
    label=r"$ m=87, p=5 $",
    color="k",
    linestyle="solid",
)
axs[2].plot(
    sr_m88p06l27[l27_index, 0],
    sr_m88p06l27[l27_index, 3],
    label=r"$ m=88, p=6 $",
    color="k",
    linestyle="dashed",
)
axs[2].plot(
    sr_m89p07l27[l27_index, 0],
    sr_m89p07l27[l27_index, 3],
    label=r"$ m=89, p=7 $",
    color="k",
    linestyle="dotted",
)

axs[1].set_xlabel(r"$N$")
axs[0].legend(loc="upper left", fontsize=7)
axs[0].set_ylabel(r"$ \epsilon_2 $")
for ax in axs:
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.hlines(e2, xmin=-60, xmax=-50, alpha=0.4)
    ax.set_xlim(-60, -50)
    ax.fill_between([-60, -50], e2_min, e2_max, alpha=0.1)

plt.savefig(fig_dir + "/e2.pdf")

# for ax in axs:
#     ax.set_ylim(1.6e-2, 2.25e-2)
#
# plt.savefig(fig_dir+"/e2_zoom.pdf")

plt.close()

fig, axs = plt.subplots(1, 3, sharey=True, sharex=True, layout="constrained")

axs[0].set_title(r"$l = 1 \times 10^{-17}$")
axs[0].plot(
    sr_m87p05l17[l17_index, 0],
    sr_m87p05l17[l17_index, 4],
    label=r"$ m=87, p=5 $",
    color="k",
    linestyle="solid",
)
axs[0].plot(
    sr_m88p06l17[l17_index, 0],
    sr_m88p06l17[l17_index, 4],
    label=r"$ m=88, p=6 $",
    color="k",
    linestyle="dashed",
)
axs[0].plot(
    sr_m89p07l17[l17_index, 0],
    sr_m89p07l17[l17_index, 4],
    label=r"$ m=89, p=7 $",
    color="k",
    linestyle="dotted",
)

axs[1].set_title(r"$l = 1 \times 10^{-22}$")
axs[1].plot(
    sr_m87p05l22[l22_index, 0],
    sr_m87p05l22[l22_index, 4],
    label=r"$ m=87, p=5 $",
    color="k",
    linestyle="solid",
)
axs[1].plot(
    sr_m88p06l22[l22_index, 0],
    sr_m88p06l22[l22_index, 4],
    label=r"$ m=88, p=6 $",
    color="k",
    linestyle="dashed",
)
axs[1].plot(
    sr_m89p07l22[l22_index, 0],
    sr_m89p07l22[l22_index, 4],
    label=r"$ m=89, p=7 $",
    color="k",
    linestyle="dotted",
)

axs[2].set_title(r"$l = 1 \times 10^{-27}$")
axs[2].plot(
    sr_m87p05l27[l27_index, 0],
    sr_m87p05l27[l27_index, 4],
    label=r"$ m=87, p=5 $",
    color="k",
    linestyle="solid",
)
axs[2].plot(
    sr_m88p06l27[l27_index, 0],
    sr_m88p06l27[l27_index, 4],
    label=r"$ m=88, p=6 $",
    color="k",
    linestyle="dashed",
)
axs[2].plot(
    sr_m89p07l27[l27_index, 0],
    sr_m89p07l27[l27_index, 4],
    label=r"$ m=89, p=7 $",
    color="k",
    linestyle="dotted",
)

axs[1].set_xlabel(r"$N$")
axs[0].set_ylabel(r"$ \epsilon_3 $")
axs[0].legend(loc="upper left", fontsize=7)
for ax in axs:
    ax.hlines(e3, xmin=-60, xmax=-50, alpha=0.4)
    ax.set_xlim(-60, -50)
    ax.fill_between([-60, -50], e3_min, e3_max, alpha=0.14)

plt.savefig(fig_dir + "/e3.pdf")

# for ax in axs:
#     ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
#     ax.set_ylim(0.015, 0.0225)
#
# plt.savefig(fig_dir+"/e3_zoom.pdf")

plt.close()
