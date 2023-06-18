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
    # lgenergy_ele = np.log10(energy_ele)  # lg(tem)
    # lgpressure_ele = np.log10(pressure_ele)  # lg(press)
    energy_ele = energy_ele.reshape(din_index, temp_index)
    pressure_ele = pressure_ele.reshape(din_index, temp_index)
    # lgenergy_ele = lgenergy_ele.reshape(din_index, temp_index)
    # lgpressure_ele = lgpressure_ele.reshape(din_index, temp_index)

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

    def get_ele_energy(ii, jj):
        e_ele = energy_ele[ii, jj] * ye
        lg_eele = math.log10(e_ele)
        return lg_eele

    def get_e_dedt(den, temp):
        din = den * ye
        ii, jj, deltalgdin, deltalgtemp = get_ii_jj(din, temp)

        minus_deltalgdin = dlgdin - deltalgdin
        minus_deltalgtemp = dlgtemp - deltalgtemp

        lgeele11 = get_ele_energy(ii, jj)
        lgeele12 = get_ele_energy(ii, jj+1)
        lgeele21 = get_ele_energy(ii+1, jj)
        lgeele22 = get_ele_energy(ii+1, jj+1)

        lgeele_t2 = deltalgdin / dlgdin * lgeele22 + minus_deltalgdin / dlgdin * lgeele12
        lgeele_t1 = deltalgdin / dlgdin * lgeele21 + minus_deltalgdin / dlgdin * lgeele11
        dlgeele_dlgt = (lgeele_t2-lgeele_t1)/dlgtemp
        lgeele_t0 = deltalgtemp / dlgtemp * lgeele_t2 + minus_deltalgtemp / dlgtemp * lgeele_t1
        e_ele = 10**lgeele_t0

        e_rad = 3 * asoli3 * temp * temp * temp * temp / den
        e_ion = 1.5 * avo * kerg * temp / a_bar
        e_total = e_ele + e_rad + e_ion

        deele_dt = e_ele*dlgeele_dlgt/temp
        deion_dt = 1.5 * avo * kerg / a_bar
        derad_dt = 12 * asoli3 * temp * temp * temp / den
        de_dt = deele_dt + deion_dt + derad_dt

        return e_total, de_dt

    def get_ele_lgpressure(ii, jj):
        p_ele = pressure_ele[ii, jj]
        lg_pele = math.log10(p_ele)
        return lg_pele

    def get_pressure(den, temp):
        din = den * ye
        ii, jj, deltalgdin, deltalgtemp = get_ii_jj(din, temp)

        minus_deltalgdin = dlgdin - deltalgdin
        minus_deltalgtemp = dlgtemp - deltalgtemp

        lgpele11 = get_ele_lgpressure(ii, jj)
        lgpele12 = get_ele_lgpressure(ii, jj+1)
        lgpele21 = get_ele_lgpressure(ii+1, jj)
        lgpele22 = get_ele_lgpressure(ii+1, jj+1)

        lgpele_t0 = deltalgtemp*minus_deltalgdin / area * lgpele12 + minus_deltalgtemp*minus_deltalgdin / area * lgpele11\
            + deltalgtemp * deltalgdin / area * lgpele22 + minus_deltalgtemp * deltalgdin / area * lgpele21
        p_ele = 10**lgpele_t0

        p_rad = asoli3 * temp * temp * temp * temp
        p_ion = den * avo * kerg * temp / a_bar
        p_total = p_ele + p_rad + p_ion

        return p_total

    count = 0
    while count < iter_times:
        e_use, dedt_use = get_e_dedt(rho, t_0)
        error_use = abs(e_use-e_want)/e_want
        if error_use < error_limit:
            break
        t_0 = t_0-(e_use-e_want)/dedt_use
        count += 1
    # print(count)

    pre_out = get_pressure(rho, t_0)

    return t_0, pre_out


print(den_ene_to_p(8.e6, 1.4e9, 2, 1, pow(10, 18), 1000, 1.e-5))
# print(den_ene_to_p(pow(10, 10.05), 1.4e9, 1, 1, pow(10, 19.05), 1000, 1.e-5))
# print(den_ene_to_p(pow(10, -2.45), 1.4e9, 1, 1, pow(10, 21.2), 1000, 1.e-5))
# print(den_ene_to_p(5.e6, 5.e9, 1, 1, 3.193786217274227e+18, 1000, 1.e-5))
