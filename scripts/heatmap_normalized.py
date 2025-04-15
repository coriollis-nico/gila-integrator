#!/usr/bin/python

# Modules
from os import makedirs
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('grayscale')
# plt.rcParams['text.usetex'] = True
# plt.rcParams['figure.figsize'] = [6.4, 14.4]

# Data import
fig_dir = "plots/heatmap"
data_dir = "data/tab/sr_variance_normalized"
makedirs(fig_dir, exist_ok=True)

# Planck data
e1_max = 0.0097
e2 = 0.032
e2_min = e2 - 0.008
e2_max = e2 + 0.009
e3 = 0.19
e3_min = e3 - 0.53
e3_max = e3 + 0.55

data = pd.read_csv(data_dir+"/variance_normalized.dat", sep=r"\s+")

data_l17 = data[(data["l"] >= 1e-18)]
data_l22 = data[(data["l"] < 1e-18) & (1e-26 < data["l"])]
data_l27 = data[(1e-26 >= data["l"])]

data_l17_mp = data_l17.loc[:, ("m", "p", "d")]
data_l17_mp_pivot = data_l17_mp.pivot(index="p",
                                      columns="m",
                                      values="d")

data_l27_mp = data_l27.loc[:, ("m", "p", "d")]
data_l27_mp_pivot = data_l27_mp.pivot(index="p",
                                      columns="m",
                                      values="d")

data_l22_mp = data_l22.loc[:, ("m", "p", "d")]
data_l22_mp_pivot = data_l22_mp.pivot(index="p",
                                      columns="m",
                                      values="d")

m_min = data_l17_mp["m"].min()
m_max = data_l17_mp["m"].max()
p_min = data_l17_mp["p"].min()
p_max = data_l17_mp["p"].max()


fig, axs = plt.subplots(1, 3, layout="constrained", squeeze=True
                        sharey=True, sharex=True)

axs[0].matshow(data_l17_mp_pivot, interpolation="none",
                    vmin=0, vmax=1,
                    extent=[m_min,
                            m_max,
                            p_max,
                            p_min])
axs[0].set_title(r"$ l = 10^{-17} $")

axs[1].matshow(data_l22_mp_pivot, interpolation="none",
                    vmin=0, vmax=1,
                    extent=[m_min,
                            m_max,
                            p_max,
                            p_min])
axs[1].set_title(r"$ l = 10^{-22} $")

axs[2].matshow(data_l27_mp_pivot, interpolation="none",
                    vmin=0, vmax=1,
                    extent=[m_min,
                            m_max,
                            p_max,
                            p_min])
axs[2].set_title(r"$ l = 10^{-27} $")


fig.colorbar(label=r"$ \sigma_1^2 \sigma_2^2 \sigma_3^2 $", extend="max")

fig.savefig(fig_dir+"/variance_normalized.pdf")
plt.close(fig)

exit()
