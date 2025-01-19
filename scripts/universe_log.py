#!/usr/bin/python

"""Finds the value of t_0 (age of the universe)
and makes a plot of a(t) using log time"""

from os import makedirs
import numpy as np
import scipy.integrate as scint
import matplotlib.pyplot as plt
plt.style.use('grayscale')
plt.rcParams['text.usetex'] = True
plt.rcParams['figure.dpi'] = 250
# plt.rcParams['figure.figsize'] = [6.4, 2.4]

fig_dir = "plots/intro_log"
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
t_tilde_i = 0
t_f = 1e-11  # /s
a_bar_i = 1
n = 200000


# Convenience functions
def rk4_step(derivative, dx, x, y):
    k1 = derivative(x, y)
    k2 = derivative(x + dx/2, y + k1 * dx/2)
    k3 = derivative(x + dx/2, y + k2 * dx/2)
    k4 = derivative(x + dx, y + k3 * dx)
    dy = dx * (k1 + 2*k2 + 2*k3 + k4)/6
    return dy


def density_bracket(a_bar, k):
    densities = O_r0/(a_bar**2) + O_m0/a_bar + k * O_k0_abs + (O_L * a_bar**2)
    return densities


tH0 = np.array(
  [scint.quad(lambda a_bar: 1/np.sqrt(
      density_bracket(a_bar, k_sign[l_])
      ), 0, 1)[0] for l_ in range(3)])
t0 = tH0 / H0


def abar_derivative(t_tilde, a_bar, l_):
    return tH0[l_] * np.exp(t_tilde) * np.sqrt(
        density_bracket(a_bar, k_sign[l_])
        )


t_tilde_f = np.log(t_f/t0)
t_tilde = np.zeros((n, 3))
dt_tilde = np.zeros(3)
for l_ in range(3):
    t_tilde[:, l_] = np.linspace(t_tilde_i, t_tilde_f[l_], n)
    dt_tilde[l_] = t_tilde[1, l_] - t_tilde[0, l_]

a_bar = np.zeros((n, 3))
a_bar[0, :] = a_bar_i

for j in range(n-1):
    for l_ in range(3):
        a_bar[j+1, l_] = a_bar[j, l_] + rk4_step(
            lambda t_tilde, a_bar: abar_derivative(t_tilde, a_bar, l_),
            dt_tilde[l_], t_tilde[j, l_], a_bar[j, l_]
            )

# a(t)
plt.figure(layout="constrained")
for l_ in range(3):
    plt.plot(np.exp(t_tilde[:, l_]), a_bar[:, l_],
             label=r"$ k = {} $".format(k_sign[l_]))
plt.xlabel(r"$ \bar{t} $")
plt.ylabel(r"$ \bar{a} $")
plt.legend(loc="best", frameon=False)
plt.savefig(fig_dir+"/scale.png", dpi=250)
plt.yscale("log", nonpositive='mask')
plt.xscale("log")
plt.savefig(fig_dir+"/scale_log.png", dpi=250)
plt.close()

# Ω_k
O_k = np.zeros((n, 3))
for l_ in range(3):
    O_k[:, l_] = k_sign[l_] * O_k0_abs/density_bracket(
        a_bar[:, l_], k_sign[l_]
        )

plt.figure(layout="constrained")
for l_ in [0, 2]:
    plt.plot(np.exp(t_tilde[:, l_]), np.abs(O_k[:, l_]),
             label=r"$ k = {} $".format(k_sign[l_]))
plt.xlabel(r"$ \bar{t} $")
plt.xscale("log")
plt.ylabel(r"$ |\Omega_{k}| $")
plt.yscale("log")
plt.legend(loc="best", frameon=False)
plt.savefig(fig_dir+"/Ok_log.png", dpi=250)
plt.close()
