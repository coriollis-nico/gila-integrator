#!/usr/bin/python

# Modules
from os import makedirs
import numpy as np
import matplotlib.pyplot as plt

plt.style.use("grayscale")
plt.rcParams["text.usetex"] = True
plt.rcParams["figure.figsize"] = [6.4, 2]
plt.rcParams["lines.linewidth"] = 1
plot_styles = ["dashed", "dashdot", "dotted"]

# Data import
fig_dir = "plots/solutions_arctanhexp"
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

data_dir = "data/sims/solutions_exp"

l_exp = np.loadtxt(data_dir + "/l.dat")

x_exp = np.loadtxt(data_dir + "/x.dat")
N_exp = np.loadtxt(data_dir + "/n.dat")

y_exp = np.loadtxt(data_dir + "/y.dat")
sr1_exp = np.loadtxt(data_dir + "/sr1.dat")
sr2_exp = np.loadtxt(data_dir + "/sr2.dat")
sr3_exp = np.loadtxt(data_dir + "/sr3.dat")

data_dir = "data/sims/solutions_arctanh"

l_ath = np.loadtxt(data_dir + "/l.dat")

x_ath = np.loadtxt(data_dir + "/x.dat")
N_ath = np.loadtxt(data_dir + "/n.dat")

y_ath = np.loadtxt(data_dir + "/y.dat")
sr1_ath = np.loadtxt(data_dir + "/sr1.dat")
sr2_ath = np.loadtxt(data_dir + "/sr2.dat")
sr3_ath = np.loadtxt(data_dir + "/sr3.dat")

# Plotting
print("Plotting…")

fig, axs = plt.subplots(nrows=1, ncols=3, squeeze=True,
                        sharex="col",
                        sharey="row",
                        layout="constrained")

for k in range(len(l_exp)):
    axs[k].set_xscale("log")
    axs[k].set_yscale("log")
    axs[k].set_title(r"$\log(l) = {}$".format(int(np.log10(l_exp[k]))))
    axs[k].plot(np.exp(x_exp),
             y_exp[k, :],
             color="k",
             ls="dashed",
             label="exp")
    axs[k].plot(np.exp(x_ath),
             y_ath[k, :],
             color="k",
             ls="dotted",
             label="atanh")

axs[1].set_xlabel(r"$ \bar{a} $")
axs[0].set_ylabel(r"$ \bar{H} $")

axs[1].legend(loc="lower left", frameon=False)

plt.savefig(fig_dir + "/solutions.png")
plt.close()


fig, axs = plt.subplots(nrows=1, ncols=3, squeeze=True,
                        sharex="col",
                        sharey="row",
                        layout="constrained")

for k in range(len(l_exp)):
    axs[k].set_xscale("log")
    axs[k].set_title(r"$\log(l) = {}$".format(int(np.log10(l_exp[k]))))
    axs[k].plot(np.exp(x_exp),
             sr1_exp[k, :],
             color="k",
             ls="dashed",
             label="exp")
    axs[k].plot(np.exp(x_ath),
             sr1_ath[k, :],
             color="k",
             ls="dotted",
             label="atanh")

axs[1].set_xlabel(r"$ \bar{a} $")
axs[0].set_ylabel(r"$ \epsilon_1 $")

axs[1].legend(loc="best", frameon=False)

plt.savefig(fig_dir + "/sr1.png")
plt.close()


fig, axs = plt.subplots(nrows=1, ncols=3, squeeze=True,
                        sharex="col",
                        sharey="row",
                        layout="constrained")

for k in range(len(l_exp)):
    axs[k].set_title(r"$\log(l) = {}$".format(int(np.log10(l_exp[k]))))
    axs[k].fill_between([N_min, N_max], 0, e1_max, alpha=0.14)
    N_index = [
        index
        for index in range(len(N_exp[k, :]))
        if (N_min <= N_exp[k, index]) and (N_exp[k, index] <= N_max)
    ]
    axs[k].plot(N_exp[k, N_index],
             sr1_exp[k, N_index],
             color="k",
             ls="dashed",
             label="exp")
    N_index = [
        index
        for index in range(len(N_ath[k, :]))
        if (N_min <= N_ath[k, index]) and (N_ath[k, index] <= N_max)
    ]
    axs[k].plot(N_ath[k, N_index],
             sr1_ath[k, N_index],
             color="k",
             ls="dotted",
             label="atanh")

axs[1].set_xlabel(r"$ \bar{a} $")
axs[0].set_ylabel(r"$ \epsilon_1 $")

axs[1].legend(loc="center right")

plt.savefig(fig_dir + "/sr1_zoom.png")
plt.close()

exit()
