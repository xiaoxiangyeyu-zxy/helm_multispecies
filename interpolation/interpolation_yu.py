import numpy as np
import math


def den_ene_to_p(rho, t_0, a_bar, z_bar, e_want, iter_times, error_limit):
    ye = z_bar/a_bar
    kerg = 1.380658e-16
    avo = 6.0221367e23
    c = 2.99792458e10
    ssol = 5.67051e-5
    asol = 4.0e0 * ssol / c
    asoli3 = asol / 3.0e0

    with open("Helmholtz_ele_106_36_den_temp.txt", "r") as yu:
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
    # lgenergy = np.log10(energy)  # lg(tem)
    # lgpressure = np.log10(pressure)  # lg(press)
    energy_ele = energy_ele.reshape(din_index, temp_index)
    pressure_ele = pressure_ele.reshape(din_index, temp_index)
    # lgenergy = lgenergy.reshape(din_index, temp_index)
    # lgpressure = lgpressure.reshape(din_index, temp_index)

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
        dlgedlgt = (lgenergy_t2-lgenergy_t1)/dlgtemp
        lgenergy_t0 = deltalgtemp / dlgtemp * lgenergy_t2 + minus_deltalgtemp / dlgtemp * lgenergy_t1
        return lgenergy_t0, dlgedlgt

    def get_total_lgpressure(ii, jj):
        din = 10**(lgdinmin+ii*dlgdin)
        den = din / ye
        temp = 10**(lgtempmin+jj*dlgtemp)
        p_ele = pressure_ele[ii, jj]
        p_rad = asoli3 * temp * temp * temp * temp
        p_ion = den * avo * kerg * temp / a_bar
        p_total = p_ele + p_rad + p_ion
        lgp_total = math.log10(p_total)
        return lgp_total

    def get_pressure(den, temp):
        din = den * ye
        ii, jj, deltalgdin, deltalgtemp = get_ii_jj(din, temp)

        minus_deltalgdin = dlgdin - deltalgdin
        minus_deltalgtemp = dlgtemp - deltalgtemp

        lgp11 = get_total_lgpressure(ii, jj)
        lgp12 = get_total_lgpressure(ii, jj+1)
        lgp21 = get_total_lgpressure(ii+1, jj)
        lgp22 = get_total_lgpressure(ii+1, jj+1)

        lgpressure_t0 = deltalgtemp*minus_deltalgdin / area * lgp12 + minus_deltalgtemp*minus_deltalgdin / area * lgp11\
            + deltalgtemp * deltalgdin / area * lgp22 + minus_deltalgtemp * deltalgdin / area * lgp21
        pressure_t0 = 10**lgpressure_t0
        return pressure_t0

    count = 0
    while count < iter_times:
        lge_use, dlgedlgt_use = get_lge_dlgedlgt(rho, t_0)
        e_use = 10**lge_use
        error_use = abs(e_use-e_want)/e_want
        if error_use < error_limit:
            break
        lg_t_0 = math.log10(t_0)-(lge_use-math.log10(e_want))/dlgedlgt_use
        t_0 = 10**lg_t_0
        count += 1
    print(count)

    pre_out = get_pressure(rho, t_0)
    return t_0, pre_out, math.log10(e_use)


print(den_ene_to_p(pow(10, 1.11), 1.4e9, 1, 1, pow(10, 19.55), 1000, 1.e-5))
print(den_ene_to_p(pow(10, 2.4), 1.4e9, 1, 1, pow(10, 20.4), 1000, 1.e-5))
print(den_ene_to_p(pow(10, 10.05), 1.4e9, 1, 1, pow(10, 19.05), 1000, 1.e-5))
print(den_ene_to_p(pow(10, -2.45), 1.4e9, 1, 1, pow(10, 21.2), 1000, 1.e-5))
print(den_ene_to_p(5.e6, 5.e9, 1, 1, 3.193786217274227e+18, 1000, 1.e-5))