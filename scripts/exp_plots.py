#!/usr/bin/python

# Modules
from os import makedirs
import numpy as np
import matplotlib.pyplot as plt

plt.style.use("grayscale")
plt.rcParams["text.usetex"] = True
# plt.rcParams["figure.figsize"] = [6.4, 2]
plt.rcParams["lines.linewidth"] = 1
plot_styles = ["dashed", "dashdot", "dotted"]

# Data import
fig_dir = "plots/solutions_exp"
data_dir = "data/sims/solutions_exp"
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

l = np.loadtxt(data_dir + "/l.dat")

x = np.loadtxt(data_dir + "/x.dat")
N = np.loadtxt(data_dir + "/n.dat")

y = np.loadtxt(data_dir + "/y.dat")
sr1 = np.loadtxt(data_dir + "/sr1.dat")
sr2 = np.loadtxt(data_dir + "/sr2.dat")
sr3 = np.loadtxt(data_dir + "/sr3.dat")

# Plotting
print("Plotting…")

plt.figure(layout="constrained")

for k in range(len(l)):
    plt.plot(np.exp(x),
             y[k, :],
             color="k",
             ls=plot_styles[k],
             label=r"$\log(l) = {}$".format(int(np.log10(l[k]))))

plt.xscale("log")
plt.yscale("log")

plt.xlim(left=np.exp(x[-1]), right=np.exp(x[0]))

plt.xlabel(r"$ \bar{a} $")
plt.ylabel(r"$ \bar{H} $")

plt.legend(loc="lower left", frameon=False)

plt.savefig(fig_dir + "/solutions.png")
plt.close()


plt.figure(layout="constrained")

for k in range(len(l)):
    plt.plot(np.exp(x),
             sr1[k, :],
             color="k",
             ls=plot_styles[k],
             label=r"$\log(l) = {}$".format(int(np.log10(l[k]))))

plt.xlim(left=np.exp(x[-1]), right=np.exp(x[0]))

plt.xscale("log")

plt.ylabel(r"$ \epsilon_1 $")
plt.xlabel(r"$ \bar{a} $")

plt.legend(loc="upper left", frameon=False)

plt.savefig(fig_dir+"/sr1.png")
plt.close()


plt.figure(layout="constrained")

for k in range(len(l)):
    N_index = [
        index
        for index in range(len(N[k, :]))
        if (N_min <= N[k, index]) and (N[k, index] <= N_max)
    ]
    plt.plot(N[k, N_index],
             sr1[k, N_index],
             color="k",
             ls=plot_styles[k],
             label=r"$\log(l) = {}$".format(int(np.log10(l[k]))))

plt.xlim(left=N_min, right=N_max)

plt.ylabel(r"$ \epsilon_1 $")
plt.xlabel(r"$ N $")

plt.fill_between([N_min, N_max], 0, e1_max, alpha=0.14)

plt.legend(loc="upper left")

plt.savefig(fig_dir+"/sr1_zoom.png")
plt.close()


plt.figure(layout="constrained")

for k in range(len(l)):
    N_index = [
        index
        for index in range(len(N[k, :]))
        if (N_min <= N[k, index]) and (N[k, index] <= N_max)
    ]
    plt.plot(N[k, N_index],
             sr2[k, N_index],
             color="k",
             ls=plot_styles[k],
             label=r"$\log(l) = {}$".format(int(np.log10(l[k]))))

plt.xlim(left=N_min, right=N_max)

plt.ylabel(r"$ \epsilon_2 $")
plt.xlabel(r"$ N $")

plt.fill_between([N_min, N_max], e2_min, e2_max, alpha=0.10)

plt.legend(loc="best")

plt.savefig(fig_dir+"/sr2.png")
plt.close()

plt.figure(layout="constrained")


for k in range(len(l)):
    N_index = [
        index
        for index in range(len(N[k, :]))
        if (N_min <= N[k, index]) and (N[k, index] <= N_max)
    ]
    plt.plot(N[k, N_index],
             sr3[k, N_index],
             color="k",
             ls=plot_styles[k],
             label=r"$\log(l) = {}$".format(int(np.log10(l[k]))))

plt.xlim(left=N_min, right=N_max)

plt.ylabel(r"$ \epsilon_3 $")
plt.xlabel(r"$ N $")

plt.fill_between([N_min, N_max], e3_min, e3_max, alpha=0.10)

plt.legend(loc="best")

plt.savefig(fig_dir+"/sr3.png")
plt.close()


exit()
