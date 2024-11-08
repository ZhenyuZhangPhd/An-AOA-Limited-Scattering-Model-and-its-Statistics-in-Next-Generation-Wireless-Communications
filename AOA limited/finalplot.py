# -*- coding: utf-8 -*-

"""
@author:zzy
@file:finalplot
@func:
@time:2022/1/17 19:42
"""
import GntData as gtd
import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':
    num = 100000
    Rr = 40
    Ri = 200
    stnum = 100
    c = 3e8
    """
    先画 CDOS模型 的 TOA 在 fig 1 中
    """
    plt.figure(1,dpi=200)
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
    plt.hist(toa, tau, density='true', color='#4196E1', histtype='step')
    # -------------------rho = Pi/40-----------------------
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
    plt.hist(toa, tau, density='true', color='#4196E1', histtype='step')
    # -------------------rho = Pi/100-----------------------
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
    plt.hist(toa, tau, density='true', color='#4196E1', histtype='step')
    # -------------------rho = Pi/0-----------------------
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
    plt.hist(toa, tau, density='true', color='#4196E1', histtype='step')
    # -------------------Draw legend-----------------------
    plt.ylim([0,1.6e7])
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    plt.legend(["rho = pi / 20", "rho = pi / 40","rho = pi / 100","rho = 0"])
    plt.title('PDF of TOA in CDOS model')
    plt.xlabel('Time of Arrival(s)')
    plt.ylabel('Probability Density')
    plt.grid(True)
    # plt.savefig("CutDiskTOAPDF.png")
    plt.show()
    print("finished")

    """
    画 CDOS模型 的 AOA 在 fig 1 中
    """
    plt.figure(3,dpi=200)

    rho = np.pi / 20
    theta_max = np.arcsin(Rr / Ri)
    theta_b = np.zeros(stnum)
    for i in range(0, stnum):
        theta_b[i] = rho + i * (theta_max - rho) / stnum
    f = gtd.CROSAOATheoreticalCurve(Ri, Rr, rho, theta_b)
    plt.plot(theta_b, f, color='k', marker='', linestyle=':')
    aoa = gtd.GetCROSAOA(Ri, 0, Rr, rho, num)
    plt.hist(aoa, theta_b, density='true', color='#4196E1', histtype='step')

    rho = np.pi / 30
    theta_b = np.zeros(stnum)
    for i in range(0, stnum):
        theta_b[i] = rho + i * (theta_max - rho) / stnum
    f2 = gtd.CROSAOATheoreticalCurve(Ri, Rr, rho, theta_b)
    plt.plot(theta_b, f2, color='k', marker='', linestyle='-.')
    aoa = gtd.GetCROSAOA(Ri, 0, Rr, rho, num)
    plt.hist(aoa, theta_b, density='true', color='#4196E1', histtype='step')

    rho = np.pi / 40
    theta_b = np.zeros(stnum)
    for i in range(0, stnum):
        theta_b[i] = rho + i * (theta_max - rho) / stnum
    f3 = gtd.CROSAOATheoreticalCurve(Ri, Rr, rho, theta_b)
    plt.plot(theta_b, f3, color='k', marker='', linestyle='--')
    aoa = gtd.GetCROSAOA(Ri, 0, Rr, rho, num)
    plt.hist(aoa, theta_b, density='true', color='#4196E1', histtype='step')

    rho = 0
    theta_b = np.zeros(stnum)
    for i in range(0, stnum):
        theta_b[i] = rho + i * (theta_max - rho) / stnum
    f4 = gtd.CROSAOATheoreticalCurve(Ri, Rr, rho, theta_b)
    plt.plot(theta_b, f4, color='k', marker='', linestyle='-')
    aoa = gtd.GetCROSAOA(Ri, 0, Rr, rho, num)
    plt.hist(aoa, theta_b, density='true', color='#4196E1', histtype='step')

    plt.legend(["rho = pi / 20", "rho = pi / 30","rho = pi / 40","rho = 0"])
    plt.title('PDF of AOA in CROS model')
    plt.xlabel('Angel of Arrival(rad)')
    plt.ylabel('Probability Density')
    plt.grid(True)
    plt.show()
    """
    画 CDOS模型 的 AOA 在 fig 1 中
    """

    plt.figure(4,dpi=200)

    rho = np.pi / 20
    theta_max = np.arcsin(Rr / Ri)
    theta_b = np.zeros(stnum)
    for i in range(0, stnum):
        theta_b[i] = rho + i * (theta_max - rho) / stnum
    f = gtd.CDOSAOATheoreticalCurve(Ri, Rr, rho, theta_b)
    plt.plot(theta_b, f, color='k', marker='', linestyle=':')
    aoa = gtd.GetCDOSAOA(Ri, 0, Rr, rho, num)
    plt.hist(aoa, theta_b, density='true', color='#4196E1', histtype='step')

    rho = np.pi / 30
    theta_b = np.zeros(stnum)
    for i in range(0, stnum):
        theta_b[i] = rho + i * (theta_max - rho) / stnum
    f2 = gtd.CDOSAOATheoreticalCurve(Ri, Rr, rho, theta_b)
    plt.plot(theta_b, f2, color='k', marker='', linestyle='-.')
    aoa = gtd.GetCDOSAOA(Ri, 0, Rr, rho, num)
    plt.hist(aoa, theta_b, density='true', color='#4196E1', histtype='step')

    rho = np.pi / 40
    theta_b = np.zeros(stnum)
    for i in range(0, stnum):
        theta_b[i] = rho + i * (theta_max - rho) / stnum
    f3 = gtd.CDOSAOATheoreticalCurve(Ri, Rr, rho, theta_b)
    plt.plot(theta_b, f3, color='k', marker='', linestyle='--')
    aoa = gtd.GetCDOSAOA(Ri, 0, Rr, rho, num)
    plt.hist(aoa, theta_b, density='true', color='#4196E1', histtype='step')

    rho = 0
    theta_b = np.zeros(stnum)
    for i in range(0, stnum):
        theta_b[i] = rho + i * (theta_max - rho) / stnum
    f4 = gtd.CDOSAOATheoreticalCurve(Ri, Rr, rho, theta_b)
    plt.plot(theta_b, f4, color='k', marker='', linestyle='-')
    aoa = gtd.GetCDOSAOA(Ri, 0, Rr, rho, num)
    plt.hist(aoa, theta_b, density='true', color='#4196E1', histtype='step')

    plt.legend(["rho = pi / 20", "rho = pi / 30", "rho = pi / 40", "rho = 0"])
    plt.title('PDF of AOA in CDOS model')
    plt.xlabel('Angel of Arrival(rad)')
    plt.ylabel('Probability Density')
    plt.grid(True)
    plt.show()






    plt.figure(2,dpi=200)
    # -------------------rho = Pi/20-----------------------
    rho = np.pi / 20
    toa = gtd.GetCROSTOA(Ri, 0, Rr, Ri, rho, num)
    tau_min = (2 * Rr + 2 * Ri * np.cos(rho) - 2 ** 0.5 * np.sqrt(
        -Ri ** 2 + 2 * Rr ** 2 + Ri ** 2 * np.cos(2 * rho))) / (2 * c)
    tau_max = (2 * Rr + 2 * Ri * np.cos(rho) + 2 ** 0.5 * np.sqrt(
        -Ri ** 2 + 2 * Rr ** 2 + Ri ** 2 * np.cos(2 * rho))) / (2 * c)
    tau = np.zeros(stnum)
    for i in range(0, stnum):
        tau[i] = tau_min + i * (tau_max - tau_min) / stnum
    f = gtd.CROSTheoreticalCurve(Ri, Rr, rho, tau)
    plt.plot(tau, f, color='k', marker='', linestyle=':')
    plt.hist(toa, tau, density='true', color='#B0C4DE', histtype='step')
    # -------------------rho = Pi/40-----------------------
    rho = np.pi / 40
    toa = gtd.GetCROSTOA(Ri, 0, Rr, Ri, rho, num)
    tau_min = (2 * Rr + 2 * Ri * np.cos(rho) - 2 ** 0.5 * np.sqrt(
        -Ri ** 2 + 2 * Rr ** 2 + Ri ** 2 * np.cos(2 * rho))) / (2 * c)
    tau_max = (2 * Rr + 2 * Ri * np.cos(rho) + 2 ** 0.5 * np.sqrt(
        -Ri ** 2 + 2 * Rr ** 2 + Ri ** 2 * np.cos(2 * rho))) / (2 * c)
    tau = np.zeros(stnum)
    for i in range(0, stnum):
        tau[i] = tau_min + i * (tau_max - tau_min) / stnum
    f2 = gtd.CROSTheoreticalCurve(Ri, Rr, rho, tau)
    plt.plot(tau, f2, color='k', marker='', linestyle='-.')
    plt.hist(toa, tau, density='true', color='#B0C4DE', histtype='step')
    # -------------------rho = Pi/100-----------------------
    rho = np.pi / 100
    toa = gtd.GetCROSTOA(Ri, 0, Rr, Ri, rho, num)
    tau_min = (2 * Rr + 2 * Ri * np.cos(rho) - 2 ** 0.5 * np.sqrt(
        -Ri ** 2 + 2 * Rr ** 2 + Ri ** 2 * np.cos(2 * rho))) / (2 * c)
    tau_max = (2 * Rr + 2 * Ri * np.cos(rho) + 2 ** 0.5 * np.sqrt(
        -Ri ** 2 + 2 * Rr ** 2 + Ri ** 2 * np.cos(2 * rho))) / (2 * c)
    tau = np.zeros(stnum)
    for i in range(0, stnum):
        tau[i] = tau_min + i * (tau_max - tau_min) / stnum
    f3 = gtd.CROSTheoreticalCurve(Ri, Rr, rho, tau)
    plt.plot(tau, f3, color='k', marker='', linestyle='--')
    plt.hist(toa, tau, density='true', color='#B0C4DE', histtype='step')
    # -------------------rho = Pi/0-----------------------
    rho = 0
    toa = gtd.GetCROSTOA(Ri, 0, Rr, Ri, rho, num)
    tau_min = (2 * Rr + 2 * Ri * np.cos(rho) - 2 ** 0.5 * np.sqrt(
        -Ri ** 2 + 2 * Rr ** 2 + Ri ** 2 * np.cos(2 * rho))) / (2 * c)
    tau_max = (2 * Rr + 2 * Ri * np.cos(rho) + 2 ** 0.5 * np.sqrt(
        -Ri ** 2 + 2 * Rr ** 2 + Ri ** 2 * np.cos(2 * rho))) / (2 * c)
    tau = np.zeros(stnum)
    for i in range(0, stnum):
        tau[i] = tau_min + i * (tau_max - tau_min) / stnum
    f4 = gtd.CROSTheoreticalCurve(Ri, Rr, rho, tau)
    plt.plot(tau, f4, color='k', marker='', linestyle='-')
    plt.hist(toa, tau, density='true', color='#B0C4DE', histtype='step')
    # -------------------Draw legend-----------------------
    plt.ylim([0, 1.6e7])
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    plt.legend(["rho = pi / 20", "rho = pi / 40","rho = pi / 100","rho = 0"])
    plt.title('PDF of TOA in CROS model')
    plt.xlabel('Time of Arrival(s)')
    plt.ylabel('Probability Density')
    plt.grid(True)
    # plt.savefig("CutDiskTOAPDF.png")
    plt.show()
    print("finished")