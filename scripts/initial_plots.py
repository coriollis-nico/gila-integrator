#!/usr/bin/python

# Modules
from os import makedirs
import numpy as np
import matplotlib.pyplot as plt

plt.style.use("seaborn-v0_8-colorblind")
plt.rcParams["text.usetex"] = True
plt.rcParams["lines.linewidth"] = 1

# Data import
fig_dir = "plots/solutions_initial"
data_dir = "data/sims/solutions_initial"
makedirs(fig_dir, exist_ok=True)

# Planck data
e1_max = 0.0097
e2 = 0.032
e2_min = e2 - 0.008
e2_max = e2 + 0.009
e3 = 0.19
e3_min = e3 - 0.53
e3_max = e3 + 0.55

# N of interest
N_max = -40
N_min = -60

# Data
print("Reading data…")

mpl = np.loadtxt(data_dir + "/mpl.dat", dtype=int)

mp = np.unique(mpl[:, 0:2], axis=0)
l = np.unique(mpl[:, 2])

x = np.loadtxt(data_dir + "/x.dat")
N = np.loadtxt(data_dir + "/n.dat")

y = np.loadtxt(data_dir + "/y.dat")
sr1 = np.loadtxt(data_dir + "/sr1.dat")
sr2 = np.loadtxt(data_dir + "/sr2.dat")
sr3 = np.loadtxt(data_dir + "/sr3.dat")
gr = np.loadtxt("data/sims/gr/gr.dat")

# Plotting
print("Plotting…")

fig, axs = plt.subplots(
    nrows=1,
    ncols=l.size,
    sharex=True,
    sharey=True,
    layout="constrained",
    squeeze=True,
    figsize=[6.4, 2.133],
)

axs[0].set_ylabel(r"$ \bar{H} $")

for k in range(l.size):
    axs[k].set_yscale("log")
    axs[k].set_xscale("log")
    axs[k].set_title(r"$ \log(l) = {} $".format(l[k]))
    axs[k].set_xlabel(r"$ \bar{a} $")
    axs[k].set_xlim(np.exp(x[-1]), np.exp(x[0]))
    axs[k].set_ylim(top=1.0e33)
    axs[k].plot(np.exp(gr[:, 0]), gr[:, 1], label="RG")

for i in range(len(mpl)):
    for k in range(l.size):
        if l[k] == mpl[i, -1]:
            this_col = k
    for j in range(len(mp)):
        if mp[j, 0] == mpl[i, 0] and mp[j, 1] == mpl[i, 1]:
            this_mp = j
    axs[this_col].plot(
        np.exp(x),
        y[i, :],
        label=r"$ m = {} $, $ p = {} $".format(mpl[i, 0], mpl[i, 1]),
    )

axs[0].legend(loc="lower left", fontsize=8, frameon=False)

plt.savefig(fig_dir + "/sols.pdf")

plt.close()

fig, axs = plt.subplots(
    nrows=3,
    ncols=l.size,
    sharex="col",
    sharey="row",
    layout="constrained",
    figsize=[6.4, 4.27],
)

for k in range(l.size):
    axs[0, k].set_title(r"$ \log(l) = {} $".format(l[k]))
    axs[2, k].set_xlabel(r"$ N $")
    for j in range(3):
        axs[j, k].set_xlim(N_min, N_max)

axs[0, 0].set_ylabel(r"$ \epsilon_1 $")
axs[1, 0].set_ylabel(r"$ \epsilon_2 $")
axs[2, 0].set_ylabel(r"$ \epsilon_3 $")

for i in range(len(mpl)):
    for k in range(l.size):
        if l[k] == mpl[i, -1]:
            this_col = k
    for j in range(len(mp)):
        if mp[j, 0] == mpl[i, 0] and mp[j, 1] == mpl[i, 1]:
            this_mp = j
    N_index = [
        index
        for index in range(len(N[i, :]))
        if (N_min <= N[i, index]) and (N[i, index] <= N_max)
    ]
    axs[0, this_col].plot(
        N[i, N_index],
        sr1[i, N_index],
        label=r"$ m = {} $, $ p = {} $".format(mpl[i, 0], mpl[i, 1]),
    )
    axs[1, this_col].plot(
        N[i, N_index],
        sr2[i, N_index],
        label=r"$ m = {} $, $ p = {} $".format(mpl[i, 0], mpl[i, 1]),
    )
    axs[2, this_col].plot(
        N[i, N_index],
        sr3[i, N_index],
        label=r"$ m = {} $, $ p = {} $".format(mpl[i, 0], mpl[i, 1]),
    )

axs[-1, 0].legend(loc="best", fontsize=6)

plt.savefig(fig_dir + "/sr_nocomp.pdf")

for k in range(l.size):
    axs[0, k].fill_between([N_min, N_max], 0, e1_max, alpha=0.2)
    axs[1, k].fill_between([N_min, N_max], e2_min, e2_max, alpha=0.1)
    axs[2, k].fill_between([N_min, N_max], e3_min, e3_max, alpha=0.2)

plt.savefig(fig_dir + "/sr.pdf")

plt.close()

exit()
