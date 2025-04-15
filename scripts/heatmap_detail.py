#!/usr/bin/python

# Modules
from os import makedirs
import pandas as pd
import matplotlib.pyplot as plt

plt.style.use("grayscale")
# plt.rcParams['text.usetex'] = True
# plt.rcParams['figure.figsize'] = [6.4, 14.4]

# Data import
fig_dir = "plots/heatmap"
data_dir = "data/tab/sr_variance_normalized"
makedirs(fig_dir, exist_ok=True)

m_min = 41
m_max = 102
p_min = 2
p_max = 4

data = pd.read_csv(data_dir + "/variance_normalized.dat", sep=r"\s+")

data_l17 = data[
    (data["l"] >= 1e-18)
    & (data["p"] <= p_max)
    & (data["p"] >= p_min)
    & (data["m"] >= m_min)
]
data_l22 = data[
    (data["l"] < 1e-18)
    & (1e-26 < data["l"])
    & (data["p"] <= p_max)
    & (data["p"] >= p_min)
    & (data["m"] >= m_min)
]
data_l27 = data[
    (1e-26 >= data["l"])
    & (data["p"] <= p_max)
    & (data["p"] >= p_min)
    & (data["m"] >= m_min)
]


data_l17_e1_mp = data_l17.loc[:, ("m", "p", "de1")]
data_l17_e1_mp_pivot = data_l17_e1_mp.pivot(index="p", columns="m", values="de1")

data_l27_e1_mp = data_l27.loc[:, ("m", "p", "de1")]
data_l27_e1_mp_pivot = data_l27_e1_mp.pivot(index="p", columns="m", values="de1")

data_l22_e1_mp = data_l22.loc[:, ("m", "p", "de1")]
data_l22_e1_mp_pivot = data_l22_e1_mp.pivot(index="p", columns="m", values="de1")


data_l17_e2_mp = data_l17.loc[:, ("m", "p", "de2")]
data_l17_e2_mp_pivot = data_l17_e2_mp.pivot(index="p", columns="m", values="de2")

data_l22_e2_mp = data_l22.loc[:, ("m", "p", "de2")]
data_l22_e2_mp_pivot = data_l22_e2_mp.pivot(index="p", columns="m", values="de2")

data_l27_e2_mp = data_l27.loc[:, ("m", "p", "de2")]
data_l27_e2_mp_pivot = data_l27_e2_mp.pivot(index="p", columns="m", values="de2")


data_l17_e3_mp = data_l17.loc[:, ("m", "p", "de3")]
data_l17_e3_mp_pivot = data_l17_e3_mp.pivot(index="p", columns="m", values="de3")

data_l22_e3_mp = data_l22.loc[:, ("m", "p", "de3")]
data_l22_e3_mp_pivot = data_l22_e3_mp.pivot(index="p", columns="m", values="de3")

data_l27_e3_mp = data_l27.loc[:, ("m", "p", "de3")]
data_l27_e3_mp_pivot = data_l27_e3_mp.pivot(index="p", columns="m", values="de3")


fig, axs = plt.subplots(3, 3, layout="constrained", sharey="col", sharex="row")

axs[0, 0].matshow(
    data_l17_e1_mp_pivot,
    interpolation="none",
    aspect="auto",
    vmin=0,
    vmax=1,
    extent=[
        m_min - 1 + 0.5,
        m_max + 0.5,
        p_max + 0.5,
        p_min - 1 + 0.5,
    ],
)
axs[0, 1].matshow(
    data_l22_e1_mp_pivot,
    interpolation="none",
    aspect="auto",
    vmin=0,
    vmax=1,
    extent=[
        m_min - 1 + 0.5,
        m_max + 0.5,
        p_max + 0.5,
        p_min - 1 + 0.5,
    ],
)
pc = axs[0, 2].matshow(
    data_l27_e1_mp_pivot,
    interpolation="none",
    aspect="auto",
    vmin=0,
    vmax=1,
    extent=[
        m_min - 1 + 0.5,
        m_max + 0.5,
        p_max + 0.5,
        p_min - 1 + 0.5,
    ],
)

fig.colorbar(pc, ax=axs[0, 2], label=r"$ \sigma_{n,1}^2 $", extend="max")

axs[1, 0].matshow(
    data_l17_e2_mp_pivot,
    interpolation="none",
    aspect="auto",
    vmin=0,
    vmax=1,
    extent=[
        m_min - 1 + 0.5,
        m_max + 0.5,
        p_max + 0.5,
        p_min - 1 + 0.5,
    ],
)
axs[1, 1].matshow(
    data_l22_e2_mp_pivot,
    interpolation="none",
    aspect="auto",
    vmin=0,
    vmax=1,
    extent=[
        m_min - 1 + 0.5,
        m_max + 0.5,
        p_max + 0.5,
        p_min - 1 + 0.5,
    ],
)
pc = axs[1, 2].matshow(
    data_l27_e2_mp_pivot,
    interpolation="none",
    aspect="auto",
    vmin=0,
    vmax=1,
    extent=[
        m_min - 1 + 0.5,
        m_max + 0.5,
        p_max + 0.5,
        p_min - 1 + 0.5,
    ],
)

fig.colorbar(pc, ax=axs[1, 2], label=r"$ \sigma_{n,2}^2 $", extend="max")

axs[2, 0].matshow(
    data_l17_e3_mp_pivot,
    interpolation="none",
    aspect="auto",
    vmin=0,
    vmax=1,
    extent=[
        m_min - 1 + 0.5,
        m_max + 0.5,
        p_max + 0.5,
        p_min - 1 + 0.5,
    ],
)
axs[2, 1].matshow(
    data_l22_e3_mp_pivot,
    interpolation="none",
    aspect="auto",
    vmin=0,
    vmax=1,
    extent=[
        m_min - 1 + 0.5,
        m_max + 0.5,
        p_max + 0.5,
        p_min - 1 + 0.5,
    ],
)
pc = axs[2, 2].matshow(
    data_l27_e3_mp_pivot,
    interpolation="none",
    aspect="auto",
    vmin=0,
    vmax=1,
    extent=[
        m_min - 1 + 0.5,
        m_max + 0.5,
        p_max + 0.5,
        p_min - 1 + 0.5,
    ],
)

fig.colorbar(pc, ax=axs[2, 2], label=r"$ \sigma_{n,3}^2 $", extend="max")


fig.savefig(fig_dir + "/variance_normalized_detail.pdf")
plt.close(fig)

exit()
