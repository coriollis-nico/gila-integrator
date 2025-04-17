#!/usr/bin/python

# Modules
from os import makedirs
import numpy as np
import matplotlib.pyplot as plt

plt.style.use("grayscale")
plt.rcParams["text.usetex"] = True
plt.rcParams["figure.figsize"] = [6.4, 2]

# Data import
fig_dir = "plots/slowroll_extra_plots"
data_dir = "data/sims/solutions"
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

mpl = np.loadtxt(data_dir + "/mpl.dat")
cases_of_interest = [0, 4, 8]

x = np.loadtxt(data_dir + "/x.dat")
N = np.loadtxt(data_dir + "/n.dat")

y = np.loadtxt(data_dir + "/y.dat")
sr1 = np.loadtxt(data_dir + "/sr1.dat")
sr2 = np.loadtxt(data_dir + "/sr2.dat")
sr3 = np.loadtxt(data_dir + "/sr3.dat")

# Plotting
print("Plotting…")

fig, axs = plt.subplots(
    nrows=1, ncols=3, sharex="col", layout="constrained", squeeze=True
)
for i in range(len(cases_of_interest)):
    N_index = [
        j for j in range(len(N[i, :])) if (N_min <= N[i, j]) and (N[i, j] <= N_max)
    ]
    axs[0].plot(
        N[i, N_index],
        sr1[i, N_index],
        label=r"$n = {n}$, $p = {p}$, $\log(l) = {ll}$".format(
            n=int(mpl[i, 0]), p=int(mpl[i, 1]), ll=int(np.log10(mpl[i, 2]))
        ),
    )
    axs[1].plot(N[i, N_index], sr2[i, N_index])
    axs[2].plot(N[i, N_index], sr3[i, N_index])

axs[0].fill_between([N_min, N_max], 0, e1_max, alpha=0.14)
axs[1].fill_between([N_min, N_max], e2_min, e2_max, alpha=0.10)
axs[2].fill_between([N_min, N_max], e3_min, e3_max, alpha=0.14)

for i in range(len(cases_of_interest)):
    axs[i].set_xlabel(r"$ N $")
    axs[i].set_xlim(N_min, N_max)

axs[0].set_ylabel(r"$ \epsilon_1 $")
axs[1].set_ylabel(r"$ \epsilon_2 $")
axs[2].set_ylabel(r"$ \epsilon_3 $")


# axs[0].legend(loc="lower left", fontsize=8)

plt.savefig(fig_dir+"/extra.pdf")

plt.close()

exit()
