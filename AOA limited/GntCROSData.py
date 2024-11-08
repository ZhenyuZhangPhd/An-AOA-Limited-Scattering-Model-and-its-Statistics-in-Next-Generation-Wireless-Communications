# -*- coding: utf-8 -*-

"""
@author:zzy
@file:GntCROSData
@func:
@time:2021/5/17 15:03
"""
import ScattererModel.GntData as gtd
import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':
    Rr = 4
    Ri = 20
    rho = np.pi / 20
    # rho = 0
    num = 100000

    # c = 3e8
    stnum = 100
    theta_max = np.arcsin(Rr / Ri)
    theta_b = np.zeros(stnum)
    for i in range(0, stnum):
        theta_b[i] = rho + i * (theta_max - rho) / stnum
    # f = gtd.CDOSAOATheoreticalCurve(Ri,Rr,rho,theta_b)
    f = gtd.CROSAOATheoreticalCurve(Ri,Rr,rho,theta_b)
    # aoa = gtd.GetCDOSAOA(Ri, 0, Rr, rho, num)
    aoa = gtd.GetCROSAOA(Ri,0,Rr,rho,num)
    # -------------------Draw curve-----------------------
    plt.plot(theta_b, f, color='k')
    plt.hist(aoa, theta_b, density='true', color='#B0C4DE')
    plt.legend(["AOA pdf", "Normalized histogram"])
    plt.title('PDF of AOA in CROS model')
    plt.xlabel('Angel of Arrival(rad)')
    plt.ylabel('Probability Density')
    plt.grid(True)
    plt.show()
