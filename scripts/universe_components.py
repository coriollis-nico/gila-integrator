#!/usr/bin/python

"Finds the value of t_0 (age of the universe) and makes a plot of a(t)"

from os import makedirs
import numpy as np
import scipy.integrate as scint
import matplotlib.pyplot as plt

plt.style.use("tableau-colorblind10")
plt.rcParams["text.usetex"] = True
plt.rcParams["figure.figsize"] = [6.4, 2.4]

fig_dir = "plots/intro"
makedirs(fig_dir, exist_ok=True)

# Integration parameters
tH0_i = 0.0
tH0_f = 1.0
n = 50000


# Functions
def a_r(tH0):
    return np.sqrt(2.0 * tH0)


def a_m(tH0):
    return (3.0 * tH0 / 2.0) ** (2 / 3)


def a_k(tH0):
    return tH0


def a_L(t):
    return 2 ** (tH0 / np.log(2)) - 1


# x axis
tH0 = np.linspace(tH0_i, tH0_f, n)

# Plotting
plt.figure(layout="constrained")
plt.xlim(0, tH0_f)
plt.plot(tH0, a_L(tH0), label=r"$ \Omega_{\Lambda} = 1 $")
plt.plot(tH0, a_m(tH0), label=r"$ \Omega_m = 1 $")
plt.plot(tH0, a_k(tH0), label=r"$ \Omega_k = 1 $")
plt.plot(tH0, a_r(tH0), label=r"$ \Omega_r = 1 $")
plt.xlabel(r"$ t H_0 $")
plt.ylabel(r"$ \bar{a}_I $")
plt.ylim(bottom=0)
plt.legend(loc="best", frameon=False)
plt.savefig(fig_dir + "/scale_component.pdf")
plt.close()
