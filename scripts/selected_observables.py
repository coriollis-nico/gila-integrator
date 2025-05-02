#!/usr/bin/python

# Modules
from os import makedirs
import numpy as np
import matplotlib.pyplot as plt

plt.style.use("grayscale")
plt.rcParams["text.usetex"] = True
plt.rcParams["lines.linewidth"] = 1
plot_styles = ["dashed", "dashdot", "dotted"]

# Data import
fig_dir = "plots/observables"
data_dir = "data/sims/solutions"
makedirs(fig_dir, exist_ok=True)

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

# Planck data
r_Pl_max = 0.063  # 95% CL
nS_Pl = 0.9668
nS_Pl_pm = 0.0037  # 68% CL
nT_Pl_min = -0.8  # 95% CL
nT_Pl_max = 0.6
alpha_S_Pl = -0.0045
alpha_S_Pl_pm = 0.0067  # 68% CL


# Observables calculations
# c.f. Arciniega et al. (2020) Inflationary predictions
#                              of geometric inflation
def r(eps_1):
    return 16 * eps_1


def nS(eps_1, eps_2, eps_3):
    C = -0.729_6
    result = (
        1
        - 2 * eps_1
        - eps_2
        - 2 * eps_1**2
        - (2 * C + 3) * eps_1 * eps_2
        - C * eps_2 * eps_3
    )
    return result


def nT(eps_1, eps_2):
    C = -0.729_6
    result = -2 * eps_1 - 2 * eps_1**2 - 2 * (C + 1) * eps_1 * eps_2
    return result


def alpha_S(eps_1, eps_2, eps_3):
    result = -2 * eps_1 * eps_2 - eps_2 * eps_3
    return result


def alpha_T(eps_1, eps_2):
    result = -2 * eps_1 * eps_2
    return result


# Plotting
print("Plotting…")

fig, axs = plt.subplots(
    nrows=5,
    ncols=l.size,
    sharex="col",
    sharey="row",
    layout="constrained",
    figsize=[6.4, 7.5],
)

for k in range(l.size):
    axs[0, k].set_title(r"$ \log(l) = {} $".format(l[k]))
    axs[2, k].set_xlabel(r"$ N $")
    for j in range(3):
        axs[j, k].set_xlim(N_min, N_max)

    axs[0, k].fill_between([N_min, N_max], 0, r_Pl_max, alpha=0.2)
    axs[1, k].fill_between(
        [N_min, N_max], nS_Pl - nS_Pl_pm, nS_Pl + nS_Pl_pm, alpha=0.1
    )
    axs[2, k].fill_between([N_min, N_max], nT_Pl_min, nT_Pl_max, alpha=0.2)
    axs[3, k].fill_between(
        [N_min, N_max],
        alpha_S_Pl - alpha_S_Pl_pm,
        alpha_S_Pl + alpha_S_Pl_pm,
        alpha=0.1,
    )

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
        r(sr1[i, N_index]),
        linestyle=plot_styles[this_mp],
        color="k",
        label=r"$ m = {} $, $ p = {} $".format(mpl[i, 0], mpl[i, 1]),
    )
    axs[1, this_col].plot(
        N[i, N_index],
        nS(sr1[i, N_index], sr2[i, N_index], sr3[i, N_index]),
        linestyle=plot_styles[this_mp],
        color="k",
        label=r"$ m = {} $, $ p = {} $".format(mpl[i, 0], mpl[i, 1]),
    )
    axs[2, this_col].plot(
        N[i, N_index],
        nT(sr1[i, N_index], sr2[i, N_index]),
        linestyle=plot_styles[this_mp],
        color="k",
        label=r"$ m = {} $, $ p = {} $".format(mpl[i, 0], mpl[i, 1]),
    )
    axs[3, this_col].plot(
        N[i, N_index],
        alpha_S(sr1[i, N_index], sr2[i, N_index], sr3[i, N_index]),
        linestyle=plot_styles[this_mp],
        color="k",
        label=r"$ m = {} $, $ p = {} $".format(mpl[i, 0], mpl[i, 1]),
    )
    axs[4, this_col].plot(
        N[i, N_index],
        alpha_T(sr1[i, N_index], sr2[i, N_index]),
        linestyle=plot_styles[this_mp],
        color="k",
        label=r"$ m = {} $, $ p = {} $".format(mpl[i, 0], mpl[i, 1]),
    )

axs[0, 0].set_ylabel(r"$ r $")
axs[1, 0].set_ylabel(r"$ n_S $")
axs[2, 0].set_ylabel(r"$ n_T $")
axs[3, 0].set_ylabel(r"$ \alpha_S $")
axs[4, 0].set_ylabel(r"$ \alpha_T $")

axs[-1, 0].legend(loc="lower left", fontsize=6, frameon=False)

plt.savefig(fig_dir + "/obs.pdf")

plt.close()

exit()
