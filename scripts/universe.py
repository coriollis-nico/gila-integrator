#!/usr/bin/python

"Finds the value of t_0 (age of the universe) and makes a plot of a(t)"

import numpy as np
import scipy.integrate as scint
import matplotlib.pyplot as plt
plt.style.use('grayscale')

# Constants
# Hubble constant / s^(-1})
H0 = 2.195e-18
# Density fractions
O_m0 = 0.3153
O_L = 0.6847
O_r0 = 9.02e-5
O_k_abs = 0.005
# (sign) of curvature
k = -1
O_k0 = k * O_k_abs
# Integration parameters
t_bar_i = 1
t_bar_f = 0
a_bar_i = 1
n = 10000
dt = (t_bar_f - t_bar_i)/n


# Convenience functions
def euler_step(derivative, dt, t, x):
    return derivative(t, x) * dt


def density_bracket(a_bar):
    densities = O_r0/(a_bar**2) + O_m0/a_bar + O_k0 + (O_L * a_bar**2)
    return np.sqrt(densities)


def seconds_to_years(t):
    return t / (1e7 * np.pi)


tH0_integral = scint.quad(lambda a_bar: 1/density_bracket(a_bar), 0, 1)
tH0 = tH0_integral[0]

t0 = tH0 / H0

print("Age of the universe:")
print(f"{seconds_to_years(t0):.3e}", "yr")
print(f"{t0:.3e}", "s")


def abar_derivative(t_bar, a_bar):
    return tH0 * density_bracket(a_bar)


a_bar = np.zeros(n)
a_bar[0] = a_bar_i

t_bar = np.zeros(n)
t_bar[0] = t_bar_i

a_bar_dot = np.zeros(n)
a_bar_dot[0] = abar_derivative(t_bar_i, a_bar_i)

for j in range(n-1):
    t_bar[j+1] = t_bar[j] + dt
    a_bar[j+1] = a_bar[j] + euler_step(abar_derivative, dt, t_bar[j], a_bar[j])
    a_bar_dot[j+1] = abar_derivative(t_bar[j+1], a_bar[j+1])

plt.figure(layout="constrained")
plt.plot(t_bar, a_bar)
plt.xlabel(r"$ \frac{t}{t_0} $")
plt.ylabel(r"$ \frac{a}{a_0} $")
plt.show()
