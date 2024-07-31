#!/usr/bin/python

"""
Plots solutions + critical point
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

fig_dir = "plots/curve_sandwich_critical"
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

sandwich_curves = [m03p01l17, m08p02l17, m10p10l17,
                   m03p01l22, m08p02l22, m10p10l22,
                   m03p01l27, m08p02l27, m10p10l27]

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

  plt.ylim(-5,75)

  plt.plot(gr_integration[0:380,0], gr_integration[0:380,1], '--', label="GR", color='0.5')

  plt.hlines( np.log(1./l[i//3]), x[0], x[-1], ls="-.", color='0.5', label=r"$ H = L^{-1} $" )

  for k in range(3):
    plt.plot(x , sandwich_curves[i+k][:,0], label="m = {}, p = {}"
             .format(mp_pairs[k,0], mp_pairs[k,1]))

    plt.legend(loc="best", frameon=False)

  plt.savefig(fig_dir+"/solcrit_{}.png".format(i//3), dpi=fig_dpi)
  plt.close()


exit()
