import numpy as np
import math

z_bar = 1
a_bar = 1
ye = z_bar/a_bar
kerg = 1.380658e-16
avo = 6.0221367e23
c = 2.99792458e10
ssol = 5.67051e-5
asol = 4.0e0 * ssol / c
asoli3 = asol / 3.0e0

with open("Helmholtz_ele_211_71_den_temp.txt", "r") as yu:
    line1 = yu.readline()
    l1 = line1.split()
    line2 = yu.readline()
    l2 = line2.split()
    data = np.loadtxt(yu)
yu.close()

din_index = int(l1[1])
dlgdin = float(l1[2])
lgdinmin = math.log10(float(l1[3]))
temp_index = int(l2[1])
dlgtemp = float(l2[2])
lgtempmin = math.log10(float(l2[3]))

area = dlgtemp*dlgdin

energy_ele = data[:, 0]  #
pressure_ele = data[:, 1]  # intensity of pressure
energy_ele = energy_ele.reshape(din_index, temp_index)
pressure_ele = pressure_ele.reshape(din_index, temp_index)


def get_ii_jj(din, temp):
    lgdin = math.log10(din)
    lgtemp = math.log10(temp)
    ii = int((lgdin-lgdinmin)//dlgdin)
    jj = int((lgtemp-lgtempmin)//dlgtemp)
    ii = max(0, min(ii, din_index-2))
    jj = max(0, min(jj, temp_index-2))
    deltalgdin = lgdin - lgdinmin - ii*dlgdin
    deltalgtemp = lgtemp - lgtempmin - jj*dlgtemp
    return ii, jj, deltalgdin, deltalgtemp


def get_total_lgenergy(ii, jj):
    din = 10**(lgdinmin+ii*dlgdin)
    den = din / ye
    temp = 10**(lgtempmin+jj*dlgtemp)
    e_ele = energy_ele[ii, jj] * ye
    e_rad = 3 * asoli3 * temp * temp * temp * temp / den
    e_ion = 1.5 * avo * kerg * temp / a_bar
    e_total = e_ele + e_rad + e_ion
    lge_total = math.log10(e_total)
    return lge_total


def get_lge_dlgedlgt(den, temp):
    din = den * ye
    ii, jj, deltalgdin, deltalgtemp = get_ii_jj(din, temp)

    minus_deltalgdin = dlgdin - deltalgdin
    minus_deltalgtemp = dlgtemp - deltalgtemp

    lge11 = get_total_lgenergy(ii, jj)
    lge12 = get_total_lgenergy(ii, jj+1)
    lge21 = get_total_lgenergy(ii+1, jj)
    lge22 = get_total_lgenergy(ii+1, jj+1)

    lgenergy_t2 = deltalgdin / dlgdin * lge22 + minus_deltalgdin / dlgdin * lge12
    lgenergy_t1 = deltalgdin / dlgdin * lge21 + minus_deltalgdin / dlgdin * lge11
    lgenergy_t0 = deltalgtemp / dlgtemp * lgenergy_t2 + minus_deltalgtemp / dlgtemp * lgenergy_t1
    return 10**lgenergy_t0


number = int(170)
den_test = 10**np.linspace(-10, 11, number)
temp_test = 10**np.linspace(4, 11, number)
e_out = [[]for i in range(number)]
for i in range(number):
    for j in range(number):
        den_use = den_test[i]
        temp_use = temp_test[j]
        e_use = get_lge_dlgedlgt(den_use, temp_use)
        e_out[i].append(e_use)

e_out = np.array(e_out)
print(e_out)
np.savetxt('C.txt', e_out)
