#!/usr/bin/python

# Modules
from os import makedirs
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use("grayscale")
plt.rcParams["text.usetex"] = True
# plt.rcParams['figure.figsize'] = [6.4, 14.4]

# Data import
fig_dir = "plots/heatmap"
data_dir = "data/tab/sr_normvar"
makedirs(fig_dir, exist_ok=True)

data = np.loadtxt(data_dir+"/vn.dat")

exit()
