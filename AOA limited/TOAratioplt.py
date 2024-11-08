# -*- coding: utf-8 -*-

"""
@author:zzy
@file:TOAratioplt
@func:
@time:2021/6/26 15:15
"""
import ScattererModel.GntData as gtd
import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':
    num = 30000
    Rr = 4
    Ri = 20
    stnum = 50
    c = 3e8
    # -------------------rho = Pi/20-----------------------
    rho = np.pi / 20
    toa = gtd.GetCDOSTOA(Ri, 0, Rr, Ri, rho, num)
    tau_min = (2 * Rr + 2 * Ri * np.cos(rho) - 2 ** 0.5 * np.sqrt(
        -Ri ** 2 + 2 * Rr ** 2 + Ri ** 2 * np.cos(2 * rho))) / (2 * c)
    tau_max = (2 * Rr + 2 * Ri * np.cos(rho) + 2 ** 0.5 * np.sqrt(
        -Ri ** 2 + 2 * Rr ** 2 + Ri ** 2 * np.cos(2 * rho))) / (2 * c)
    tau = np.zeros(stnum)
    for i in range(0, stnum):
        tau[i] = tau_min + i * (tau_max - tau_min) / stnum
    f = gtd.CDOSTheoreticalCurve(Ri, Rr, rho, tau)
    plt.plot(tau, f, color='k', marker='', linestyle=':')
    # -------------------rho2-----------------------
    rho = np.pi / 40
    toa = gtd.GetCDOSTOA(Ri, 0, Rr, Ri, rho, num)
    tau_min = (2 * Rr + 2 * Ri * np.cos(rho) - 2 ** 0.5 * np.sqrt(
        -Ri ** 2 + 2 * Rr ** 2 + Ri ** 2 * np.cos(2 * rho))) / (2 * c)
    tau_max = (2 * Rr + 2 * Ri * np.cos(rho) + 2 ** 0.5 * np.sqrt(
        -Ri ** 2 + 2 * Rr ** 2 + Ri ** 2 * np.cos(2 * rho))) / (2 * c)
    tau = np.zeros(stnum)
    for i in range(0, stnum):
        tau[i] = tau_min + i * (tau_max - tau_min) / stnum
    f2 = gtd.CDOSTheoreticalCurve(Ri, Rr, rho, tau)
    plt.plot(tau, f2, color='k', marker='', linestyle='-.')
    # -------------------Draw curve-----------------------
    rho = np.pi / 100
    toa = gtd.GetCDOSTOA(Ri, 0, Rr, Ri, rho, num)
    tau_min = (2 * Rr + 2 * Ri * np.cos(rho) - 2 ** 0.5 * np.sqrt(
        -Ri ** 2 + 2 * Rr ** 2 + Ri ** 2 * np.cos(2 * rho))) / (2 * c)
    tau_max = (2 * Rr + 2 * Ri * np.cos(rho) + 2 ** 0.5 * np.sqrt(
        -Ri ** 2 + 2 * Rr ** 2 + Ri ** 2 * np.cos(2 * rho))) / (2 * c)
    tau = np.zeros(stnum)
    for i in range(0, stnum):
        tau[i] = tau_min + i * (tau_max - tau_min) / stnum
    f3 = gtd.CDOSTheoreticalCurve(Ri, Rr, rho, tau)
    plt.plot(tau, f3, color='k', marker='', linestyle='--')

    rho = 0
    toa = gtd.GetCDOSTOA(Ri, 0, Rr, Ri, rho, num)
    tau_min = (2 * Rr + 2 * Ri * np.cos(rho) - 2 ** 0.5 * np.sqrt(
        -Ri ** 2 + 2 * Rr ** 2 + Ri ** 2 * np.cos(2 * rho))) / (2 * c)
    tau_max = (2 * Rr + 2 * Ri * np.cos(rho) + 2 ** 0.5 * np.sqrt(
        -Ri ** 2 + 2 * Rr ** 2 + Ri ** 2 * np.cos(2 * rho))) / (2 * c)
    tau = np.zeros(stnum)
    for i in range(0, stnum):
        tau[i] = tau_min + i * (tau_max - tau_min) / stnum
    f4 = gtd.CDOSTheoreticalCurve(Ri, Rr, rho, tau)
    plt.plot(tau, f4, color='k', marker='', linestyle='-')
    # -------------------Draw curve-----------------------
    plt.legend(["rho = pi / 20", "rho = pi / 40","rho = pi / 100","rho = 0"])
    plt.title('PDF of TOA in CDOS model')
    plt.xlabel('Time of Arrival(s)')
    plt.ylabel('Probability Density')
    plt.grid(True)
    # plt.savefig("CutDiskTOAPDF.png")
    plt.show()
    print("finished")