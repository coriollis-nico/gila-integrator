#!/usr/bin/python

# Modules
from os import makedirs
import numpy as np
import matplotlib.pyplot as plt

plt.style.use("grayscale")
plt.rcParams["text.usetex"] = True
plt.rcParams["figure.figsize"] = [6.4, 6.4]
import pandas as pd

# Data import
fig_dir = "plots/heatmap"
makedirs(fig_dir, exist_ok=True)
data_dir = "data/tab/sr_normvar"

# Create tab
data = pd.read_csv(data_dir + "/vn.dat", sep=r"\s+")

# Parameters
m = np.unique(data["m"].to_numpy())
p = np.unique(data["p"].to_numpy())
l = np.unique(data["log_l"].to_numpy())


# Plot
fig, axs = plt.subplots(
    nrows=3, ncols=l.size, sharex="all", sharey="all", layout="constrained"
)

v = ["v1", "v2", "v3"]
for k in range(l.size):
    axs[0, k].set_title(r"$ \log l = {} $".format(l[k]))
    axs[2, k].set_xlabel(r"$ m $")
    for vc in range(len(v)):
        axs[vc, k].matshow(
            data[data["log_l"] == l[k]].pivot_table(
                values=v[vc], index="p", columns="m"
            ),
            vmin=0,
            vmax=1,
            interpolation="none",
            extent=(m[0] - 0.5, m[-1] - 0.5, p[-1] - 0.5, p[0] - 0.5),
            aspect="auto",
        )

        axs[k, 0].set_ylabel(r"$ p $")

plt.savefig(fig_dir + "/heatmap.pdf")

plt.close()

exit()
