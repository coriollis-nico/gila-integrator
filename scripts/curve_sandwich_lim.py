#!/usr/bin/python

"""
Plots solutions + critical point + limit curves
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

fig_dir = "plots/curve_sandwich_lim"
data_dir = "data/sims/curve_sandwich"
makedirs(fig_dir, exist_ok = True)

## 0)

notification("Reading `gr.dat`")
gr_integration = np.loadtxt("data/sims/gr/gr.dat", comments='#')

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

l17 = np.loadtxt(data_dir+"/lim_l1.0E-17.dat", comments='#')
l22 = np.loadtxt(data_dir+"/lim_l1.0E-22.dat", comments='#')
l27 = np.loadtxt(data_dir+"/lim_l1.0E-27.dat", comments='#')

sandwich_curves = [m03p01l17, m08p02l17, m10p10l17,
                   m03p01l22, m08p02l22, m10p10l22,
                   m03p01l27, m08p02l27, m10p10l27]
limit_curves = [l17, l22, l27]

e1_limits = (-1, 0.97e-2) # actually lower is -inf
e2_limits = (2.6e-2, 4.1e-2)
e3_limits = (-340e-3, 640e-3)

mp_pairs = np.loadtxt(data_dir+"/mp.dat", comments="#", dtype=int)
l = np.loadtxt(data_dir+"/l.dat", comments="#")

notification("Creating plots...")

for i in range(0,9,3):

  # Solutions
  plt.figure()

  plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

  plt.title(r"Solutions (l = {})".format(l[i//3]))

  plt.xlabel(r"$ \ln{\frac{a}{a_0}} $")
  plt.ylabel(r"$ \ln{\frac{H}{H_0}} $")

  plt.ylim(-1,70)

  plt.plot(gr_integration[0:380,0], gr_integration[0:380,1], ':', label="GR", color='0.75')

  plt.plot(x, limit_curves[i//3][:,0], label=r"$ m,p \to \infty $", color='k', lw=0.9)

  for k in range(3):
    plt.plot(x , sandwich_curves[i+k][:,0], ls='--', label="m = {}, p = {}"
             .format(mp_pairs[k,0], mp_pairs[k,1]))

  plt.legend(loc="best", frameon=False)

  plt.savefig(fig_dir+"/solcrit_{}.svg".format(i//3))
  plt.close()

  # slow-roll 1
  plt.figure()

  plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

  plt.title("Slow-roll 1 (l = {})".format(l[i//3]))

  plt.xlabel(r"$ \ln{\frac{a}{a_0}} $")
  plt.ylabel(r"$ \epsilon_1 $")

  plt.ylim(-0.2e-2,0.037)

  plt.plot(x[zoom], limit_curves[i//3][zoom,1], label=r"$ m,p \to \infty $", color='k', lw=0.9)

  for k in range(3):
    plt.plot(x[zoom] , sandwich_curves[i+k][zoom,1], ls='--', label="m = {}, p = {}"
             .format(mp_pairs[k,0], mp_pairs[k,1]))

  plt.legend(loc="best", frameon=False)

  plt.fill_between(x[zoom], e1_limits[0], e1_limits[1], label="Planck constraint (95% CL)",
                   alpha = 0.15, color = 'purple', lw = 0)

  plt.savefig(fig_dir+"/e1crit_{}.svg".format(i//3))
  plt.close()

  ## unzoom
  plt.figure()

  plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

  plt.title("Slow-roll 1 (l = {})".format(l[i//3]))

  plt.xlabel(r"$ \ln{\frac{a}{a_0}} $")
  plt.ylabel(r"$ \epsilon_1 $")

  plt.plot(x, limit_curves[i//3][:,1], label=r"$ m,p \to \infty $", color='k', lw=0.9)

  for k in range(3):
    plt.plot(x , sandwich_curves[i+k][:,1], ls='--', label="m = {}, p = {}"
             .format(mp_pairs[k,0], mp_pairs[k,1]))

  plt.legend(loc="best", frameon=False)

  plt.savefig(fig_dir+"/e1uzcrit_{}.svg".format(i//3))
  plt.close()

  # slow-roll 2
  plt.figure()

  plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

  plt.title("Slow-roll 2 (l = {})".format(l[i//3]))

  plt.xlabel(r"$ \ln{\frac{a}{a_0}} $")
  plt.ylabel(r"$ \epsilon_2 $")

  plt.ylim(2.6e-2, 7e-2)

  for k in range(3):
    plt.plot(x[zoom] , sandwich_curves[i+k][zoom,2], ls='--', label="m = {}, p = {}"
             .format(mp_pairs[k,0], mp_pairs[k,1]))

  plt.legend(loc="best", frameon=False)

  plt.fill_between(x[zoom], e2_limits[0], e2_limits[1], label="Planck constraint (95% CL)",
                   alpha = 0.15, color = 'purple', lw = 0)

  plt.savefig(fig_dir+"/e2crit_{}.svg".format(i//3))
  plt.close()

  # slow-roll 3
  plt.figure()

  plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

  plt.title("Slow-roll 3 (l = {})".format(l[i//3]))

  plt.xlabel(r"$ \ln{\frac{a}{a_0}} $")
  plt.ylabel(r"$ \epsilon_3 $")

  plt.ylim(2.5e-2,6.5e-2)

  for k in range(3):
    plt.plot(x[zoom] , sandwich_curves[i+k][zoom,3], ls='--', label="m = {}, p = {}"
             .format(mp_pairs[k,0], mp_pairs[k,1]))

  plt.legend(loc="best", frameon=False)

  plt.fill_between(x[zoom], e3_limits[0], e3_limits[1], label="Planck constraint (95% CL)",
                   alpha = 0.15, color = 'purple', lw = 0)

  plt.savefig(fig_dir+"/e3crit_{}.svg".format(i//3))
  plt.close()

exit()
