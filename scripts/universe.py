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
t_bar_i = 1
t_bar_f = 0
a_bar_i = 1
n = 50000
dt = (t_bar_f - t_bar_i) / n


# Convenience functions
def euler_step(derivative, dt, t, x):
    return derivative(t, x) * dt


def density_bracket(a_bar, k):
    densities = O_r0 / (a_bar**2) + O_m0 / a_bar + k * O_k0_abs + (O_L * a_bar**2)
    return np.sqrt(densities)


def seconds_to_years(t):
    return t / (1e7 * np.pi)


tH0 = np.array(
    [scint.quad(lambda a_bar: 1 / density_bracket(a_bar, k), 0, 1)[0] for k in k_sign]
)
t0 = tH0 / H0

print("Age of the universe:")
for j in range(len(k_sign)):
    print(
        k_sign[j], "->", f"{t0[j]:.3e}", "s =", f"{seconds_to_years(t0[j]):.3e}", "yr"
    )


def abar_derivative(t, a, i):
    return tH0[i] * density_bracket(a, k_sign[i])


t_bar = np.zeros(n)
t_bar[0] = t_bar_i

a_bar = np.zeros((n, len(k_sign)))
a_bar[0, :] = a_bar_i

d_abar_d_tbar = np.zeros((n, len(k_sign)))
for l_ in range(len(k_sign)):
    d_abar_d_tbar[0, l_] = abar_derivative(t_bar_i, a_bar_i, l_)

for j in range(n - 1):
    t_bar[j + 1] = t_bar[j] + dt
    for l_ in range(len(k_sign)):
        a_bar[j + 1, l_] = a_bar[j, l_] + euler_step(
            (lambda t, a: abar_derivative(t, a, k_sign[l_])), dt, t_bar[j], a_bar[j, l_]
        )
        d_abar_d_tbar[j + 1, l_] = abar_derivative(t_bar[j + 1], a_bar[j + 1, l_], l_)

# a(t)
plt.figure(layout="constrained")
plt.xlim(0, 1)
for l_ in range(len(k_sign)):
    plt.plot(t_bar, a_bar[:, l_], label=r"$ k = {} $".format(k_sign[l_]))
plt.xlabel(r"$ \bar{t} $")
plt.ylabel(r"$ \bar{a} $")
plt.legend(loc="best", frameon=False)
plt.savefig(fig_dir + "/scale.pdf")
plt.close()

# 1/a_dot
plt.figure(layout="constrained")
plt.xlim(0, 1)
for l_ in range(len(k_sign)):
    plt.plot(t_bar, 1 / d_abar_d_tbar[:, l_], label=r"$ k = {} $".format(k_sign[l_]))
plt.xlabel(r"$ \bar{t} $")
plt.ylabel(r"$ \left( \frac{\mathrm{d} \bar{a}}{\mathrm{d} \bar{t}} \right)^{-1} $")
plt.legend(loc="best", frameon=False)
plt.savefig(fig_dir + "/scale_deriv_inverse.pdf")
plt.close()

# Δa
plt.figure(layout="constrained")
plt.xlim(0, 1)
for l_ in [0, 2]:
    plt.plot(
        t_bar,
        (a_bar[:, l_] - a_bar[:, 1]) / a_bar[:, 1],
        label=r"$ k = {} $".format(k_sign[l_]),
    )
plt.xlabel(r"$ \bar{t} $")
plt.ylabel(r"$ \frac{\Delta{\bar{a}}}{\bar{a}_{k=0}} $")
plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
plt.legend(loc="best", frameon=False)
plt.savefig(fig_dir + "/scale_diff.pdf")
plt.close()

# Ω_k
O_k = np.zeros((n, len(k_sign)))
for l_ in range(len(k_sign)):
    O_k[:, l_] = k_sign[l_] * O_k0_abs * (tH0[l_] / d_abar_d_tbar[:, l_]) ** 2

print("Last values of Ω_k:")
print("t =", f"{t_bar[-1]:.3e}")
for l_ in [0, 2]:
    print(k_sign[l_], "->", f"{O_k[-1, l_]:.3e}")

plt.figure(layout="constrained")
plt.xlim(0, 1)
for l_ in [0, 2]:
    plt.plot(t_bar, O_k[:, l_], label=r"$ k = {} $".format(k_sign[l_]))
plt.xlabel(r"$ \bar{t} $")
plt.ylabel(r"$ \Omega_{k} $")
plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
plt.legend(loc="best", frameon=False)
plt.savefig(fig_dir + "/O_k.pdf")
plt.close()
