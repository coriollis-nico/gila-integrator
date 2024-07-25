#!/usr/bin/python

"""
- Changes sign of curve_sandwich slow rolls & compares with Planck Data
"""

import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn-v0_8-colorblind')

from os import makedirs

quiet = False

def notification(info):
  if (quiet == False):
    print(info)

fig_dpi = 300

fig_dir = "plots/switched_intervals"
data_dir = "data/sims/curve_sandwich"
makedirs(fig_dir, exist_ok = True)

## 0)

notification("Reading sandwich_curves")
x = np.loadtxt(data_dir+"/x_grid.dat", comments='#')
zoom = ( x > -60 ) & ( x < -50 )

# Columnas: y // e_j
#           0     j(->4)

m03p01l17 = np.loadtxt(data_dir+"/m03p01l1.0E-17.dat", comments='#')
m03p01l22 = np.loadtxt(data_dir+"/m03p01l1.0E-22.dat", comments='#')
m03p01l27 = np.loadtxt(data_dir+"/m03p01l1.0E-27.dat", comments='#')

m08p02l17 = np.loadtxt(data_dir+"/m08p02l1.0E-17.dat", comments='#')
m08p02l22 = np.loadtxt(data_dir+"/m08p02l1.0E-22.dat", comments='#')
m08p02l27 = np.loadtxt(data_dir+"/m08p02l1.0E-27.dat", comments='#')

m10p10l17 = np.loadtxt(data_dir+"/m10p10l1.0E-17.dat", comments='#')
m10p10l22 = np.loadtxt(data_dir+"/m10p10l1.0E-22.dat", comments='#')
m10p10l27 = np.loadtxt(data_dir+"/m10p10l1.0E-27.dat", comments='#')

sandwich_curves = [m03p01l17, m08p02l17, m10p10l17,
                   m03p01l22, m08p02l22, m10p10l22,
                   m03p01l27, m08p02l27, m10p10l27]

mp_pairs = np.loadtxt(data_dir+"/mp.dat", comments="#", dtype=int)
l = np.loadtxt(data_dir+"/l.dat", comments="#")


notification("Creating plots...")

for i in range(0,9,3):

  # slow-roll 1
  plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

  plt.title("1st slow-roll parameter (switched) (l = {})".format(l[i//3]))

  plt.ylim(0, 3.6e-2)

  plt.xlabel(r"$ \ln{\frac{a}{a_0}} $")
  plt.ylabel(r"$ \epsilon_1 $")

  for k in range(3):
    plt.plot(x[zoom], -sandwich_curves[i+k][zoom,1], label="m = {}, p = {}"
             .format(mp_pairs[k,0], mp_pairs[k,1]))

    plt.legend(loc="best")

  plt.fill_between(x[zoom], 9.7e-3, 0, label="Planck constraint (95% CL)",
                   alpha = 0.15, color = 'purple', lw = 0)

  plt.savefig(fig_dir+"/e1zs_{}.png".format(i//3), dpi=fig_dpi)
  plt.close()

  # slow-roll 2
  plt.figure()

  plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

  plt.title("2st slow-roll parameter (switched) (l = {})".format(l[i//3]))

  plt.ylim(2.5e-2, 6.9e-2)

  plt.xlabel(r"$ \ln{\frac{a}{a_0}} $")
  plt.ylabel(r"$ \epsilon_2 $")

  for k in range(3):
    plt.plot(x[zoom], -sandwich_curves[i+k][zoom,2], label="m = {}, p = {}"
             .format(mp_pairs[k,0], mp_pairs[k,1]))

    plt.legend(loc="best")

  plt.fill_between(x[zoom], 4.1e-2, 2.6e-2, label="Planck constraint (68% CL)",
                   alpha = 0.15, color = 'red', lw = 0)

  plt.savefig(fig_dir+"/e2zs_{}.png".format(i//3), dpi=fig_dpi)
  plt.close()

  # slow-roll 3
  plt.figure()

  plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

  plt.title("3rd slow-roll parameter (switched) (l = {})".format(l[i//3]))

  plt.ylim(0, 4.8e-3)

  plt.xlabel(r"$ \ln{\frac{a}{a_0}} $")
  plt.ylabel(r"$ \epsilon_3 $")

  for k in range(3):
    plt.plot(x[zoom], -sandwich_curves[i+k][zoom,3], label="m = {}, p = {}"
             .format(mp_pairs[k,0], mp_pairs[k,1]))

    plt.legend(loc="best")

  plt.fill_between(x[zoom], 640e-3, -340e-3, label="Planck constraint (95% CL)",
                   alpha = 0.15, color = 'purple', lw = 0)

  plt.savefig(fig_dir+"/e3zs_{}.png".format(i//3), dpi=fig_dpi)
  plt.close()

exit()
