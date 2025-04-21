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
fig_dir = "plots/observables"
data_dir = "data/sims/solutions"
makedirs(fig_dir, exist_ok=True)

# N of interest
N_max = -40
N_min = -60

# Data
print("Reading data…")

mpl = np.loadtxt(data_dir + "/mpl.dat")

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


def alpha_T(eps_1, eps_2, eps_3):
    result = -2 * eps_1 * eps_2
    return result


# Plotting
print("Plotting…")

fig, axs = plt.subplots(
    nrows=5,
    ncols=3,
    sharex="col",
    sharey="row",
    layout="constrained",
    figsize=[6.4, 7.5],
)

for k in range(3):
    axs[0, k].set_title(r"$\log(l) = {}$".format(int(np.log10(mpl[k, 2]))))

    axs[4, k].set_xlabel(r"$ N $")

    axs[0, k].fill_between([N_min, N_max], 0, r_Pl_max, alpha=0.14)
    axs[1, k].fill_between(
        [N_min, N_max], nS_Pl - nS_Pl_pm, nS_Pl + nS_Pl_pm, alpha=0.10
    )
    axs[2, k].fill_between([N_min, N_max], nT_Pl_min, nT_Pl_max, alpha=0.14)
    axs[3, k].fill_between(
        [N_min, N_max],
        alpha_S_Pl - alpha_S_Pl_pm,
        alpha_S_Pl + alpha_S_Pl_pm,
        alpha=0.10,
    )

    for j in range(4):
        axs[j, k].set_xlim(N_min, N_max)

    for i in range(3):
        N_index = [
            l
            for l in range(len(N[3 * i + k, :]))
            if (N_min <= N[3 * i + k, l]) and (N[3 * i + k, l] <= N_max)
        ]

        axs[0, k].plot(
            N[3 * i + k, N_index],
            r(sr1[3 * i + k, N_index]),
            color="k",
            ls=plot_styles[i],
            label=r"$n = {}$, $p = {}$".format(
                int(mpl[3 * i + k, 0]), int(mpl[3 * i + k, 1])
            ),
        )
        axs[1, k].plot(
            N[3 * i + k, N_index],
            nS(
                sr1[3 * i + k, N_index],
                sr2[3 * i + k, N_index],
                sr3[3 * i + k, N_index],
            ),
            color="k",
            ls=plot_styles[i],
            label=r"$n = {}$, $p = {}$".format(
                int(mpl[3 * i + k, 0]), int(mpl[3 * i + k, 1])
            ),
        )
        axs[2, k].plot(
            N[3 * i + k, N_index],
            nT(sr1[3 * i + k, N_index], sr2[3 * i + k, N_index]),
            color="k",
            ls=plot_styles[i],
            label=r"$n = {}$, $p = {}$".format(
                int(mpl[3 * i + k, 0]), int(mpl[3 * i + k, 1])
            ),
        )
        axs[3, k].plot(
            N[3 * i + k, N_index],
            alpha_S(
                sr1[3 * i + k, N_index],
                sr2[3 * i + k, N_index],
                sr3[3 * i + k, N_index],
            ),
            color="k",
            ls=plot_styles[i],
            label=r"$n = {}$, $p = {}$".format(
                int(mpl[3 * i + k, 0]), int(mpl[3 * i + k, 1])
            ),
        )
        axs[4, k].plot(
            N[3 * i + k, N_index],
            alpha_T(
                sr1[3 * i + k, N_index],
                sr2[3 * i + k, N_index],
                sr3[3 * i + k, N_index],
            ),
            color="k",
            ls=plot_styles[i],
            label=r"$n = {}$, $p = {}$".format(
                int(mpl[3 * i + k, 0]), int(mpl[3 * i + k, 1])
            ),
        )

axs[0, 0].set_ylabel(r"$ r $")
axs[1, 0].set_ylabel(r"$ n_S $")
axs[2, 0].set_ylabel(r"$ n_T $")
axs[3, 0].set_ylabel(r"$ \alpha_S $")
axs[4, 0].set_ylabel(r"$ \alpha_T $")

axs[2, 1].legend(loc="best", fontsize=8)

plt.savefig(fig_dir + "/obs.pdf")

plt.close()

exit()
