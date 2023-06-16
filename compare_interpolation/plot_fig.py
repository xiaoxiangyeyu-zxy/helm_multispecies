import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, colors


def plot_fig(ee, x1, x2, name):
    plt.figure()
    xp1, yp1 = x1, x2
    plt.contourf(xp1, yp1, ee)
    plt.xlabel('lg_density')
    plt.ylabel('lg_temperature')
    plt.title('error_'+name)
    plt.colorbar()
    plt.savefig(name + '.png')

