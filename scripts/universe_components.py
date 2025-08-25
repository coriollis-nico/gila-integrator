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

# Integration parameters
t_i = 0.0
t_f = 1.0
n = 50000

# Functions
def a_r(t):
    return np.sqrt(t)

def a_m(t):
    return t**(2/3)

def a_k(t):
    return t

def a_L(t):
    return 2**(t) - 1

# x axis
t = np.linspace(t_i, t_f, n)

# Plotting
plt.figure(layout="constrained")

plt.plot(t, a_L(t), color="k", ls="solid", label=r"$ \Omega_{\Lambda} = 1 $")
plt.plot(t, a_m(t), color="k", ls="dashed", label=r"$ \Omega_m = 1 $")
plt.plot(t, a_k(t), color="k", ls="dotted", label=r"$ \Omega_k = 1 $")
plt.plot(t, a_r(t), color="k", ls="dashdot", label=r"$ \Omega_r = 1 $")

plt.xlabel(r"$ \bar{t} $")
plt.ylabel(r"$ \bar{a}_I $")
plt.legend(loc="best", frameon=False)
plt.savefig(fig_dir + "/scale_component.pdf")
plt.close()
