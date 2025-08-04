#!/usr/bin/python

"Finds the value of t_0 (age of the universe) and makes a plot of a(t)"

from os import makedirs
import numpy as np
import scipy.integrate as scint
import matplotlib.pyplot as plt

plt.style.use("grayscale")
plt.rcParams["text.usetex"] = True
plt.rcParams["figure.figsize"] = [6.4, 2.4]

fig_dir = "plots/intro"
makedirs(fig_dir, exist_ok=True)

# Constants
# Hubble constant / s^(-1})
H0 = 2.195e-18
# Density fractions
O_m0 = 0.3153
O_L = 0.6847
O_r0 = 9.02e-5
O_k0_abs = 0.005
k_sign = np.array([-1, 0, 1])
# Integration parameters
t_bar_i = 1.0
t_bar_f = 0.0
a_bar_i = 1.0
n = 50000


# Convenience functions
def seconds_to_years(t):
    return t / (1e7 * np.pi)


# t0H0
tH0 = np.array(
    [
        scint.quad(
            lambda a: 1.0 / np.sqrt(O_r0 / a**2 + O_m0 / a + k * O_k0_abs + O_L * a**2),
            0,
            1,
        )[0]
        for k in k_sign
    ]
)
t0 = tH0 / H0

print("Age of the universe:")
for j in range(len(k_sign)):
    print(
        k_sign[j], "->", f"{t0[j]:.3e}", "s =", f"{seconds_to_years(t0[j]):.3e}", "yr"
    )


# deriavatives
def da(t, a, i):
    return tH0[i] * np.sqrt(O_r0 / a**2 + O_m0 / a + k_sign[i] * O_k0_abs + O_L * a**2)


def da2(t, a, i):
    return -(0.5) * (tH0[i] ** 2) * (2.0 * O_r0 / a**3 + O_m0 / a**2 - 2.0 * O_L * a)


# data arrays
t_bar = np.linspace(t_bar_i, t_bar_f, n)
dt = t_bar[1] - t_bar[0]

a_bar = np.zeros((n, len(k_sign)))
a_bar[0, :] = a_bar_i

adot = np.zeros((n, len(k_sign)))
for k_i in range(len(k_sign)):
    adot[0, k_i] = da(t_bar[0], a_bar[0, k_i], k_i)

addot = np.zeros((n, len(k_sign)))
for k_i in range(len(k_sign)):
    addot[0, k_i] = da2(t_bar[0], a_bar[0, k_i], k_i)


# integration
for k_i in range(len(k_sign)):
    for n_i in range(n - 1):
        a_bar[n_i + 1, k_i] = a_bar[n_i, k_i] + (dt * adot[n_i, k_i])
        adot[n_i + 1, k_i] = da(t_bar[n_i + 1], a_bar[n_i + 1, k_i], k_i)
        addot[n_i + 1, k_i] = da2(t_bar[n_i + 1], a_bar[n_i + 1, k_i], k_i)


# a(t)
plt.figure(layout="constrained")
plt.xlim(0, 1)
for k_i in range(len(k_sign)):
    plt.plot(t_bar, a_bar[:, k_i], label=r"$ k = {} $".format(k_sign[k_i]))
plt.xlabel(r"$ \bar{t} $")
plt.ylabel(r"$ \bar{a} $")
plt.legend(loc="best", frameon=False)
plt.savefig(fig_dir + "/scale.pdf")
plt.close()

# d2a
plt.figure(layout="constrained")
for k_i in range(len(k_sign)):
    plt.plot(t_bar, addot[:, k_i], label=r"$ k = {} $".format(k_sign[k_i]))
plt.legend(loc="best", frameon=False)
plt.xlim(0, 0.1)
plt.yscale("symlog")
plt.ylim(top=0)
plt.xlabel(r"$ \bar{t} $")
plt.ylabel(r"$ \ddot{\bar{a}} $")
plt.savefig(fig_dir + "/ddot_scale.pdf")
plt.close()


# Δa
# plt.figure(layout="constrained")
# plt.xlim(0, 1)
# for k_i in range(len(k_sign)):
#     plt.plot(
#         t_bar,
#         (a_bar[:, k_i] - a_bar[:, 1]) / a_bar[:, 1],
#         label=r"$ k = {} $".format(k_sign[k_i]),
#     )
# plt.xlabel(r"$ \bar{t} $")
# plt.ylabel(r"$ \frac{\Delta{\bar{a}}}{\bar{a}_{k=0}} $")
# plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
# plt.legend(loc="best", frameon=False)
# plt.savefig(fig_dir + "/scale_diff.pdf")
# plt.close()
