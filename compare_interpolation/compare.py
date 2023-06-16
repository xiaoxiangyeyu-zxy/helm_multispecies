import numpy as np
import math
from plot_fig import plot_fig

A = np.loadtxt("A.txt")
B = np.loadtxt("B.txt")
C = np.loadtxt("C.txt")
i_5 = np.loadtxt("5.txt")

error_A = abs((A-i_5)/i_5)
error_B = abs((B-i_5)/i_5)
error_C = abs((C-i_5)/i_5)
errormean_A = np.mean(error_A)
errormean_B = np.mean(error_B)
errormean_C = np.mean(error_C)
print(errormean_A)
print(errormean_B)
print(errormean_C)

number = 170
dd = 21/(number-1)
dt = 7/(number-1)
x_1, x_2 = np.mgrid[-10: 11+dd: dd, 4: 11+dt: dt]

plot_fig(error_A, x_1, x_2, 'A')
plot_fig(error_B, x_1, x_2, 'B')
plot_fig(error_C, x_1, x_2, 'C')
