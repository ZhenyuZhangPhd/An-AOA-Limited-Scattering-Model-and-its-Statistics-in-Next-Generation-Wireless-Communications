# -*- coding: utf-8 -*-

"""
@author:zzy
@file:ROSdata
@func:该主程序生成TOA的直方图和理论曲线，可以调整为其他模型
@time:2021/5/12 20:28

"""
import numpy as np
import CutDiskTOAPDF as cdt
import CutDiskAOAPDF as cda
import matplotlib.pyplot as plt


def GetAOD(xy, Ri):
    """
    xy表示散射体坐标，在该函数中，假设基站坐标为(0,0)，假设移动台坐标为(Ri,0)
    :param xy:散射体坐标
    :param Ri:基站离移动台的实际距离
    :return:返回该由散射体得到的TOA
    """
    l = np.size(xy[1])
    aod = np.zeros(l)
    for i in range(0, l):
        aod[i] = np.arctan(xy[1, i] / (Ri - xy[0, i]))
    return aod


def GetAOA(xy):
    """
    xy表示散射体坐标，在该函数中，假设基站坐标为(0,0)，假设移动台坐标为(Ri,0)
    :param xy:散射体坐标
    :return:返回该由散射体得到的TOA
    """
    l = np.size(xy[1])
    aoa = np.zeros(l)
    for i in range(0, l):
        aoa[i] = np.arctan(xy[1, i] / xy[0, i])
    return aoa


def GetROSTOA(Rx, Ry, Radius, Ri, num):
    """

    :param Rx:一般取值为Ri,表示移动台离基站的x距离
    :param Ry:一般取值为0,表示移动台离基站的y距离
    :param Radius:散射体半径的平方，一般取值为Rr
    :param Ri: 基站离移动台的实际距离
    :param num: 仿真生成独立同分布的样本数量
    :return: 返回该由散射体得到的TOA
    """
    xy = np.zeros((2, num))
    for i in range(0, num):
        xy[:, i] = cdt.UniCutXY(Rx, Ry, Radius, rho=0, Disk_or_Ring=0)
    toa = cdt.GetTOA(xy, Ri)
    return toa


def GetDOSTOA(Rx, Ry, Radius, Ri, num):
    xy = np.zeros((2, num))
    for i in range(0, num):
        xy[:, i] = cdt.UniCutXY(Rx, Ry, Radius, rho=0, Disk_or_Ring=1)
    toa = cdt.GetTOA(xy, Ri)
    return toa


def GetCDOSAOA(Rx, Ry, Radius, rho, num):
    xy = np.zeros((2, num))
    for i in range(0, num):
        xy[:, i] = cdt.UniCutXY(Rx, Ry, Radius, rho, 1)
    aoa = cda.GetAOA(xy)
    return aoa


def GetCROSAOA(Rx, Ry, Radius, rho, num):
    xy = np.zeros((2, num))
    for i in range(0, num):
        xy[:, i] = cdt.UniCutXY(Rx, Ry, Radius, rho, Disk_or_Ring=0)
    aoa = cda.GetAOA(xy)
    return aoa


def GetCDOSTOA(Rx, Ry, Radius, Ri, rho, num):
    """

    :param Rx:一般取值为Ri,表示移动台离基站的x距离
    :param Ry:一般取值为0,表示移动台离基站的y距离
    :param Radius:散射体半径的平方，一般取值为Rr
    :param Ri: 基站离移动台的实际距离
    :param rho:CDOS模型中的遮挡角度大小，用弧度表示，例如np.pi / 20
    :param num: 仿真生成独立同分布的样本数量
    :return: 返回该由散射体得到的TOA
    """
    xy = np.zeros((2, num))
    for i in range(0, num):
        xy[:, i] = cdt.UniCutXY(Rx, Ry, Radius, rho, 1)
    toa = cdt.GetTOA(xy, Ri)
    return toa


def GetCROSTOA(Rx, Ry, Radius, Ri, rho, num):
    xy = np.zeros((2, num))
    for i in range(0, num):
        xy[:, i] = cdt.UniCutXY(Rx, Ry, Radius, rho, Disk_or_Ring=0)
    toa = cdt.GetTOA(xy, Ri)
    return toa


def CDOSTheoreticalCurve(Ri, Rr, rho, tau):
    c = 3e8
    # A = (Rr ** 2) * (2 * np.arccos(np.sin(rho) * Ri / Rr) - np.sin(2 * np.arccos(np.sin(rho) * Ri / Rr)))
    # A这里应该没有问题，是化简后的结果
    A = (Rr ** 2) * (np.pi-2 * np.arcsin(np.sin(rho) * Ri / Rr) - np.sin(2 * np.arcsin(np.sin(rho) * Ri / Rr)))
    f1 = np.sin(rho) * (c * tau * np.cos(rho) - Ri) * Ri ** 2 / (c * tau - np.cos(rho) * Ri) ** 2
    f2 = np.pi * (Ri ** 2 - 2 * (c * tau) ** 2) / np.sqrt((c * tau) ** 2 - Ri ** 2)
    f3 = 2 * np.arctan(np.sqrt((c * tau) ** 2 - Ri ** 2) * np.tan(rho / 2) / (Ri - c * tau)) * (
            Ri ** 2 - 2 * (c * tau) ** 2) / np.sqrt((c * tau) ** 2 - Ri ** 2)
    f4 = 2 * np.arctan(np.sqrt((c * tau) ** 2 - Ri ** 2) * np.tan(
        0.5 * np.arccos((Ri ** 2 + c * tau * (2 * Rr - tau * c)) / (2 * Ri * Rr))) / (Ri - c * tau)) * (
                 Ri ** 2 - 2 * c ** 2 * tau ** 2) / np.sqrt(c ** 2 * tau ** 2 - Ri ** 2)
    f5 = Ri * Rr * (2 * Rr - c * tau) * np.sqrt(
        (c ** 2 * tau ** 2 - Ri ** 2) * (Ri ** 2 - (c * tau - 2 * Rr) ** 2) / (Ri ** 2 * Rr ** 2)) / (
                 Ri ** 2 - c ** 2 * tau ** 2)
    f = (c / (4 * A)) * (f1 - f2 - f3 - f4 - f5)
    return f


def CDOSAOATheoreticalCurve(Ri, Rr, rho, theta_b):
    A = (Rr ** 2) * (2 * np.arccos(np.sin(rho) * Ri / Rr) - np.sin(2 * np.arccos(np.sin(rho) * Ri / Rr)))
    f = 2 * Ri * np.cos(theta_b) * np.sqrt((Ri ** 2) * (np.cos(theta_b) ** 2) - Ri ** 2 + Rr ** 2) / A
    f = f * 2
    return f


def CROSTheoreticalCurve(Ri, Rr, rho, tau):
    c = 3e8
    f1=c*(tau*c-Rr)
    f2=((Ri**2 + Rr**2 -(tau*c-Rr)**2)**2)/(4* Ri**2 * Rr**2)
    f3=Ri*Rr*np.sqrt(1-f2)

    f4=np.pi-2*np.arcsin(Ri*np.sin(rho)/Rr)

    f = (f1/f3)/f4
    return f
# def ROSTheroratialCurve(Ri,Rr,tau):


def CROSAOATheoreticalCurve(Ri, Rr, rho,theta_b):
    X = np.pi-2*np.arcsin(Ri*np.sin(rho)/Rr)
    A = np.sqrt(-(Ri ** 2) / 2 + Rr ** 2 + 0.5 * Ri ** 2 * np.cos(2 * theta_b))
    B = np.sqrt(Rr ** 2 + Ri ** 2 * np.cos(2 * theta_b) - Ri * np.cos(theta_b) * 2 * A)
    C = np.sqrt(Rr ** 2 + Ri ** 2 * np.cos(2 * theta_b) + Ri * np.cos(theta_b) * 2 * A)
    f = (B+C)/(X*A)
    return f

if __name__ == '__main__':
    Rr = 4
    Ri = 20
    rho = np.pi / 20
    # rho=0
    num = 10000
    # toa = np.zeros((1,num))
    # toa = GetCDOSTOA(Ri, 0, Rr, Ri, rho, num)
    toa = GetCROSTOA(Ri,0,Rr,Ri,rho,num)
    # toa = GetROSTOA(Ri,0,Rr,Ri,num)
    # toa = GetDOSTOA(Ri,0,Rr,Ri,num)
    c = 3e8
    stnum = 50
    tau_min = (2 * Rr + 2 * Ri * np.cos(rho) - 2 ** 0.5 * np.sqrt(
        -Ri ** 2 + 2 * Rr ** 2 + Ri ** 2 * np.cos(2 * rho))) / (2 * c)
    tau_max = (2 * Rr + 2 * Ri * np.cos(rho) + 2 ** 0.5 * np.sqrt(
        -Ri ** 2 + 2 * Rr ** 2 + Ri ** 2 * np.cos(2 * rho))) / (2 * c)
    tau = np.zeros(stnum)

    for i in range(0, stnum):
        tau[i] = tau_min + i * (tau_max - tau_min) / stnum

    f = CROSTheoreticalCurve(Ri, Rr, rho, tau)
    # -------------------Draw curve-----------------------
    plt.figure(figsize=(10, 5))

    plt.subplot(1,2,1)
    plt.plot(tau, f, color='k')
    plt.hist(toa, tau, density='true', color='#B0C4DE',histtype='step')
    plt.legend(["TOA pdf", "Normalized histogram"])
    plt.title('PDF of TOA')
    plt.xlabel('Time of Arrival(s)')
    plt.ylabel('Probability Density')
    plt.grid(False)

    plt.subplot(1, 2, 2)
    plt.plot(tau, f, color='k')
    plt.hist(toa, tau, density='true', color='#B0C4DE')
    plt.legend(["TOA pdf", "Normalized histogram"])
    plt.title('PDF of TOA')
    plt.xlabel('Time of Arrival(s)')
    plt.ylabel('Probability Density')
    plt.grid(True)
    # plt.savefig("CutDiskTOAPDF.png")
    plt.show()
    print("finished")
