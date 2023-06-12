import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, colors


def plot_fig(ee, pp, x1, x2, x_shape, y_shape):
    lev_exp_ee = np.linspace(np.floor(np.log10(ee.min())), np.ceil(np.log10(ee.max())), 40)
    levs_ee = np.power(10, lev_exp_ee)

    plt.figure()
    xp1, yp1 = x1, x2
    zp1 = ee.reshape(x_shape, y_shape)
    plt.contourf(xp1, yp1, zp1, levs_ee, norm=colors.LogNorm())
    plt.xlabel('density')
    plt.ylabel('temperature')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('energy_ele')
    plt.colorbar(format=ticker.LogFormatter(10, labelOnlyBase=False))
    plt.savefig('energy.png')

    lev_exp_pp = np.linspace(np.floor(np.log10(pp.min())), np.ceil(np.log10(pp.max())), 40)
    levs_pp = np.power(10, lev_exp_pp)

    plt.figure()
    xp2, yp2 = x1, x2
    zp2 = pp.reshape(x_shape, y_shape)
    plt.contourf(xp2, yp2, zp2, levs_pp, norm=colors.LogNorm())
    plt.xlabel('density')
    plt.ylabel('temperature')
    plt.xscale('log')
    plt.yscale('log')
    plt.title('pressure_ele')
    plt.colorbar(format=ticker.LogFormatter(10, labelOnlyBase=False))
    plt.savefig('pressure.png')
