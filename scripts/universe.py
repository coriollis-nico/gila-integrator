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
k_sign = np.array([1, 0, -1])
# Integration parameters
t_i = 1.0
t_l = 0.0
t_r = 2
a_i = 1.0
n_i = 50000
n = 2 * n_i


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
        k_sign[j], "->", f"{t0[j]:.3e}", "m =", f"{seconds_to_years(t0[j]):.3e}", "yr"
    )


# deriavatives
def da(t, a, i):
    return tH0[i] * np.sqrt(O_r0 / a**2 + O_m0 / a + k_sign[i] * O_k0_abs + O_L * a**2)


def da2(t, a, i):
    return -(0.5) * (tH0[i] ** 2) * (2.0 * O_r0 / a**3 + O_m0 / a**2 - 2.0 * O_L * a)


# data arrays
t = np.linspace(t_l, t_r, n)
dt = t[1] - t[0]

a = np.zeros((n, len(k_sign)))
a[n_i, :] = a_i

adot = np.zeros((n, len(k_sign)))
for k_i in range(len(k_sign)):
    adot[n_i, k_i] = da(t[n_i], a[n_i, k_i], k_i)

addot = np.zeros((n, len(k_sign)))
for k_i in range(len(k_sign)):
    addot[n_i, k_i] = da2(t[n_i], a[n_i, k_i], k_i)


# integration
for k_i in range(len(k_sign)):
    # future
    for j in range(1, n - n_i):
        a[n_i + j, k_i] = a[n_i + j - 1, k_i] + (dt * adot[n_i + j - 1, k_i])
        adot[n_i + j, k_i] = da(t[n_i + j], a[n_i + j, k_i], k_i)
        addot[n_i + j, k_i] = da2(t[n_i + j], a[n_i + j, k_i], k_i)
    # past
    for j in range(1, n_i + 1):
        a[n_i - j, k_i] = a[n_i - j + 1, k_i] - (dt * adot[n_i - j + 1, k_i])
        adot[n_i - j, k_i] = da(t[n_i - j], a[n_i - j, k_i], k_i)
        addot[n_i - j, k_i] = da2(t[n_i - j], a[n_i - j, k_i], k_i)

# a(t)
plt.figure(layout="constrained")
plt.xlim(0, 2)
for k_i in range(len(k_sign)):
    plt.plot(t, a[:, k_i], label=r"$ k = {} $".format(k_sign[k_i]))
plt.xlabel(r"$ \bar{t} $")
plt.ylabel(r"$ \bar{a} $")
# plt.yscale("log")
plt.legend(loc="best", frameon=False)
plt.savefig(fig_dir + "/scale.pdf")
plt.close()

# relative densities
a_X = np.logspace(-5, 1, n)

rho_m_norm = O_m0 / (a_X**3)
rho_r_norm = O_r0 / (a_X**4)
rho_L_norm = O_L / (a_X**0)
rho_k_norm_abs = O_k0_abs / (a_X**2)

plt.figure(layout="constrained")
plt.xscale("log")
plt.yscale("log")
plt.plot(a_X, rho_L_norm, label=r"$ \rho_{\Lambda} $")
plt.plot(a_X, rho_m_norm, label=r"$ \rho_{m} $")
plt.plot(a_X, rho_k_norm_abs, label=r"$ \rho_{k} $ ($ k = -1 $)")
plt.plot(a_X, rho_r_norm, label=r"$ \rho_{r} $")
plt.xlim(1e-5, 10)
plt.xlabel(r"$ \bar{a} $")
plt.ylabel(r"$ \frac{\rho_I}{\rho_{c,0}} $")
plt.legend(loc="best", frameon=False)
plt.savefig(fig_dir + "/contribution.pdf")
plt.close()


# d2a
plt.figure(layout="constrained")
for k_i in range(len(k_sign)):
    plt.plot(t, addot[:, k_i], label=r"$ k = {} $".format(k_sign[k_i]))
plt.legend(loc="best", frameon=False)
plt.xlim(0, 0.1)
plt.yscale("symlog")
plt.ylim(top=0)
plt.xlabel(r"$ \bar{t} $")
plt.ylabel(r"$ \ddot{\bar{a}} $")
plt.savefig(fig_dir + "/ddot_scale.pdf")
plt.close()
