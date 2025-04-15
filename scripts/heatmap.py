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
data_dir = "data/tab/sr_variance"
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

data_l17 = data[(data["l"] >= 1e-18)]
data_l22 = data[(data["l"] < 1e-18) & (1e-26 < data["l"])]
data_l27 = data[(1e-26 >= data["l"])]


data_l17_e1_mp = data_l17.loc[:, ("m", "p", "Δϵ1")]
data_l17_e1_mp_pivot = data_l17_e1_mp.pivot(index="p",
                                      columns="m",
                                      values="Δϵ1")

data_l27_e1_mp = data_l27.loc[:, ("m", "p", "Δϵ1")]
data_l27_e1_mp_pivot = data_l27_e1_mp.pivot(index="p",
                                      columns="m",
                                      values="Δϵ1")

data_l22_e1_mp = data_l22.loc[:, ("m", "p", "Δϵ1")]
data_l22_e1_mp_pivot = data_l22_e1_mp.pivot(index="p",
                                      columns="m",
                                      values="Δϵ1")


data_l17_e2_mp = data_l17.loc[:, ("m", "p", "σ²2")]
data_l17_e2_mp_pivot = data_l17_e2_mp.pivot(index="p",
                                      columns="m",
                                      values="σ²2")

data_l22_e2_mp = data_l22.loc[:, ("m", "p", "σ²2")]
data_l22_e2_mp_pivot = data_l22_e2_mp.pivot(index="p",
                                      columns="m",
                                      values="σ²2")

data_l27_e2_mp = data_l27.loc[:, ("m", "p", "σ²2")]
data_l27_e2_mp_pivot = data_l27_e2_mp.pivot(index="p",
                                      columns="m",
                                      values="σ²2")


data_l17_e3_mp = data_l17.loc[:, ("m", "p", "σ²3")]
data_l17_e3_mp_pivot = data_l17_e3_mp.pivot(index="p",
                                      columns="m",
                                      values="σ²3")

data_l22_e3_mp = data_l22.loc[:, ("m", "p", "σ²3")]
data_l22_e3_mp_pivot = data_l22_e3_mp.pivot(index="p",
                                      columns="m",
                                      values="σ²3")

data_l27_e3_mp = data_l27.loc[:, ("m", "p", "σ²3")]
data_l27_e3_mp_pivot = data_l27_e3_mp.pivot(index="p",
                                      columns="m",
                                      values="σ²3")


fig, axs = plt.subplots(3, 3, layout="constrained",
                        sharey="col", sharex="row")

axs[0, 0].matshow(data_l17_e1_mp_pivot, interpolation="none",
                    vmin=-e1_max, vmax=0,
                    extent=[data_l17_e1_mp["m"].min(),
                            data_l17_e1_mp["m"].max(),
                            data_l17_e1_mp["p"].max(),
                            data_l17_e1_mp["p"].min()])
axs[0, 1].matshow(data_l22_e1_mp_pivot, interpolation="none",
                    vmin=-e1_max, vmax=0,
                    extent=[data_l22_e1_mp["m"].min(),
                            data_l22_e1_mp["m"].max(),
                            data_l22_e1_mp["p"].max(),
                            data_l22_e1_mp["p"].min()])
axs[0, 2].matshow(data_l27_e1_mp_pivot, interpolation="none",
                    vmin=-e1_max, vmax=0,
                    extent=[data_l27_e1_mp["m"].min(),
                            data_l27_e1_mp["m"].max(),
                            data_l27_e1_mp["p"].max(),
                            data_l27_e1_mp["p"].min()])

axs[1, 0].matshow(data_l17_e2_mp_pivot, interpolation="none",
                    vmin=0, vmax=0.009**2,
                    extent=[data_l17_e2_mp["m"].min(),
                            data_l17_e2_mp["m"].max(),
                            data_l17_e2_mp["p"].max(),
                            data_l17_e2_mp["p"].min()])
axs[1, 1].matshow(data_l22_e2_mp_pivot, interpolation="none",
                    vmin=0, vmax=0.009**2,
                    extent=[data_l22_e2_mp["m"].min(),
                            data_l22_e2_mp["m"].max(),
                            data_l22_e2_mp["p"].max(),
                            data_l22_e2_mp["p"].min()])
axs[1, 2].matshow(data_l27_e2_mp_pivot, interpolation="none",
                    vmin=0, vmax=0.009**2,
                    extent=[data_l27_e2_mp["m"].min(),
                            data_l27_e2_mp["m"].max(),
                            data_l27_e2_mp["p"].max(),
                            data_l27_e2_mp["p"].min()])

axs[2, 0].matshow(data_l17_e3_mp_pivot, interpolation="none",
                    vmin=0, vmax=0.55**2,
                    extent=[data_l17_e3_mp["m"].min(),
                            data_l17_e3_mp["m"].max(),
                            data_l17_e3_mp["p"].max(),
                            data_l17_e3_mp["p"].min()])
axs[2, 1].matshow(data_l22_e3_mp_pivot, interpolation="none",
                    vmin=0, vmax=0.55**2,
                    extent=[data_l22_e3_mp["m"].min(),
                            data_l22_e3_mp["m"].max(),
                            data_l22_e3_mp["p"].max(),
                            data_l22_e3_mp["p"].min()])
axs[2, 2].matshow(data_l27_e3_mp_pivot, interpolation="none",
                    vmin=0, vmax=0.55**2,
                    extent=[data_l27_e3_mp["m"].min(),
                            data_l27_e3_mp["m"].max(),
                            data_l27_e3_mp["p"].max(),
                            data_l27_e3_mp["p"].min()])


# fig.colorbar(pc, label=r"$ \Delta \epsilon_1 $", extend="max", shrink=0.5)

fig.savefig(fig_dir+"/variance.pdf")
plt.close(fig)

exit()
