# -*- coding: utf-8 -*-

"""
@author:zzy
@file:ratioplt
@func:
@time:2021/6/25 10:27
"""
import ScattererModel.GntData as gtd
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    num = 30000
    Rr = 4
    Ri = 20

    rho = np.pi / 20
    # rho = 0
    stnum = 50
    theta_max = np.arcsin(Rr / Ri)
    theta_b = np.zeros(stnum)
    for i in range(0, stnum):
        theta_b[i] = rho + i * (theta_max - rho) / stnum
    f = gtd.CROSAOATheoreticalCurve(Ri,Rr,rho,theta_b)
    plt.plot(theta_b, f, color='k', marker='', linestyle=':')

    rho = np.pi / 30
    theta_b = np.zeros(stnum)
    for i in range(0, stnum):
        theta_b[i] = rho + i * (theta_max - rho) / stnum
    f2 = gtd.CROSAOATheoreticalCurve(Ri, Rr, rho, theta_b)
    plt.plot(theta_b, f2, color='k', marker='', linestyle='-.')

    rho = np.pi / 40
    theta_b = np.zeros(stnum)
    for i in range(0, stnum):
        theta_b[i] = rho + i * (theta_max - rho) / stnum
    f3 = gtd.CROSAOATheoreticalCurve(Ri, Rr, rho, theta_b)
    plt.plot(theta_b, f3, color='k', marker='', linestyle='--')

    rho = 0
    theta_b = np.zeros(stnum)
    for i in range(0, stnum):
        theta_b[i] = rho + i * (theta_max - rho) / stnum
    f4 = gtd.CROSAOATheoreticalCurve(Ri, Rr, rho, theta_b)
    plt.plot(theta_b, f4, color='k', marker='', linestyle='-')

    plt.legend(["rho = pi / 20", "rho = pi / 30","rho = pi / 40","rho = 0"])
    plt.title('PDF of AOA in CROS model')
    plt.xlabel('Angel of Arrival(rad)')
    plt.ylabel('Probability Density')
    plt.grid(True)
    plt.show()