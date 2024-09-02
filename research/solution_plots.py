#!/usr/bin/python

"""
Plots solutions
"""

# Modules
from os import makedirs
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn-v0_8-colorblind')

# Data import
fig_dir = "plots/solution_plots"
data_dir = "data/sims/curve_sandwich"
makedirs(fig_dir, exist_ok=True)


print("Reading `gr.dat`")
gr_integration = np.loadtxt("data/sims/gr/gr.dat", comments='#')

print("Reading sandwich_curves")
x = np.loadtxt(data_dir+"/x_grid.dat", comments='#')

m03p01l17 = np.loadtxt(data_dir+"/m03p01l1.0E-17.dat", comments='#')
m03p01l22 = np.loadtxt(data_dir+"/m03p01l1.0E-22.dat", comments='#')
m03p01l27 = np.loadtxt(data_dir+"/m03p01l1.0E-27.dat", comments='#')

m08p02l17 = np.loadtxt(data_dir+"/m08p02l1.0E-17.dat", comments='#')
m08p02l22 = np.loadtxt(data_dir+"/m08p02l1.0E-22.dat", comments='#')
m08p02l27 = np.loadtxt(data_dir+"/m08p02l1.0E-27.dat", comments='#')

m10p10l17 = np.loadtxt(data_dir+"/m10p10l1.0E-17.dat", comments='#')
m10p10l22 = np.loadtxt(data_dir+"/m10p10l1.0E-22.dat", comments='#')
m10p10l27 = np.loadtxt(data_dir+"/m10p10l1.0E-27.dat", comments='#')

l17 = np.loadtxt(data_dir+"/lim_l1.0E-17.dat", comments='#')
l22 = np.loadtxt(data_dir+"/lim_l1.0E-22.dat", comments='#')
l27 = np.loadtxt(data_dir+"/lim_l1.0E-27.dat", comments='#')

sandwich_curves = [m03p01l17, m08p02l17, m10p10l17,
                   m03p01l22, m08p02l22, m10p10l22,
                   m03p01l27, m08p02l27, m10p10l27]
limit_curves = [l17, l22, l27]

mp_pairs = np.loadtxt(data_dir+"/mp.dat", comments="#", dtype=int)

l_value = np.loadtxt(data_dir+"/l.dat", comments="#")

# Plotting

plt.figure()
plt.yscale("log")
plt.xscale("log")
plt.xlabel(r"$ \frac{a}{a_0} $")
plt.ylabel(r"$ \frac{H}{H_0} $")
plt.plot(np.exp(x), sandwich_curves[0])
plt.show()

exit()
