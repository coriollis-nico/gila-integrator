#!/usr/bin/python

"""
Plots the values of l as a function of when the inflation exit is expected to happen
"""

import matplotlib.pyplot as plt
import numpy as np
from os import makedirs
plt.style.use('grayscale')

fig_dir = "plots/inflation_exit"
makedirs(fig_dir, exist_ok = True)

# function
def l(a, Omega_M, Omega_R, Omega_Dark):
  return 1./np.sqrt(
    Omega_M/2./a**3 + Omega_R/2./a**4 + (3./2.)*Omega_Dark*np.log(a)
    - (Omega_M + Omega_R)/2. + 1
    )

# a_bar values
a_bar = np.logspace(-44,0)
## Density values
Omega_M0 = 0.31
Omega_R0 = 8.4e-5
Omega_D0 = 0.69

plt.figure()
plt.title("Salida de inflación")
plt.ticklabel_format(style='sci', scilimits=(0,0))
plt.xscale("log")
plt.yscale("log")

plt.xlabel(r"$ \bar{a}_I $")
plt.ylabel(r"$ l^{-1} $")

plt.plot( a_bar, 1./l(a_bar, Omega_M0, Omega_R0, Omega_D0))
plt.show()
