#!/usr/bin/python

# Modules
from os import makedirs
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('grayscale')
# plt.rcParams['text.usetex'] = True
plt.rcParams['figure.figsize'] = [6.4, 3]

# Data import
fig_dir = "plots/sorting"
data_dir = "data/tab/more_cases"
makedirs(fig_dir, exist_ok=True)

# Planck data
e1_max = 0.0097
e2 = 0.032
e2_min = e2 - 0.008
e2_max = e2 + 0.009
e3 = 0.19
e3_min = e3 - 0.53
e3_max = e3 + 0.55

data = pd.read_csv(data_dir+"/variance.dat", sep=r"\s+")

data_l22 = data[(data["l"] < 1e-18) & (1e-26 < data["l"])]

data_l22_mp = data_l22.loc[:, ("m", "p", "σ²2")]
data_l22_mp_pivot = data_l22_mp.pivot(index="p",
                                      columns="m",
                                      values="σ²2")
data_l22_mp_matrix = data_l22_mp_pivot.to_numpy()

# NOTE: Asume que no hay hoyos en lista de m o de p

print(data_l22_mp_pivot)
print(data_l22_mp_matrix)

plt.figure(layout="constrained")

plt.matshow(data_l22_mp_matrix, fignum=0,
            vmax=0.009**2,
            extent=[data_l22_mp["m"].min(),
                    data_l22_mp["m"].max(),
                    data_l22_mp["p"].max(),
                    data_l22_mp["p"].min()])

plt.colorbar()

plt.show()

# data_l22_mp.plot.scatter(x="m", y="p", c="σ²2",
#                          marker="s", s=190,
#                          cmap="grey",
#                          vmax=9e-3,  # vmin=np.nanmin(data_l22_mp_matrix),
#                          norm="log")
#
# plt.xlim(left=2.5, right=52.5)
# plt.ylim(0.5, 50.5)
#
# plt.show()
#
# plt.close()

exit()
