import numpy as np
import math
from plot_fig import plot_fig

# read the data from the file "Helmholtz_ele_211_71_den_temp.txt"
with open("Helmholtz_ele_211_71_den_temp.txt", "r") as f:
    line1 = f.readline()
    l1 = line1.split()
    line2 = f.readline()
    l2 = line2.split()
    data = np.loadtxt(f)
f.close()

den_index = int(l1[1])
dDen = float(l1[2])
lgdenmin = math.log10(float(l1[3]))
temp_index = int(l2[1])
dT = float(l2[2])
lgtempmin = math.log10(float(l2[3]))

print(data)

# the first column of data is temperature, the second column of data is intensity of pressure
E = data[:, 0]
P = data[:, 1]

# the first index is density, the second index is energy
Euse = E.reshape(den_index, temp_index)
Puse = P.reshape(den_index, temp_index)

lgdenmax = lgdenmin + dDen*(den_index-1)
lgtempmax = lgtempmin + dT*(temp_index-1)
x_shape = den_index
y_shape = temp_index

x1, x2 = np.mgrid[lgdenmin:(lgdenmax+dDen):dDen, lgtempmin:(lgtempmax+dT):dT]  # 141*101
x1_plot = 10**x1
x2_plot = 10**x2

plot_fig(Euse, Puse, x1_plot, x2_plot, x_shape, y_shape)
