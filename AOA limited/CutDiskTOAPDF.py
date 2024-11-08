import numpy as np
import matplotlib.pyplot as plt


def UniCutXY(Rx, Ry, Radius, rho, Disk_or_Ring):
    """

    :param Disk_or_Ring:若是1 表示圆盘，0表示圆环
    :param Rx:一般取值为Ri,表示移动台离基站的x距离
    :param Ry:一般取值为0,表示移动台离基站的y距离
    :param Radius:散射体半径的平方，一般取值为Rr//没有平方
    :param rho:CDOS模型中的遮挡角度大小，用弧度表示，例如np.pi / 20
    :return:返回生成散射体的xy坐标，是ndarray(2,1)的向量，其中xy[0, 0]表示x轴，xy[1, i]表示y轴
    """
    while (1):
        theta = np.random.random() * 2 * np.pi
        R = 0
        if Disk_or_Ring == 1:
            R = np.random.uniform(0, Radius ** 2)
            R = R ** 0.5
        elif Disk_or_Ring == 0:
            R = Radius
        x = np.sin(theta) * R
        y = np.cos(theta) * R
        if (((y + Ry) >= ((x + Rx) * np.tan(rho))) | ((y + Ry) <= ((x + Rx) * np.tan(-rho)))):
            u = np.array((x + Rx, y + Ry))
            break
    return u


def GetTOA(xy, Ri):
    """
    xy表示散射体坐标，在该函数中，假设基站坐标为(0,0)，假设移动台坐标为(Ri,0)
    :param xy:散射体坐标
    :param Ri:基站离移动台的实际距离
    :return:返回该由散射体得到的TOA
    """
    l = np.size(xy[1])
    toa = np.zeros(l)
    for i in range(0, l):
        toa[i] = (np.sqrt(xy[0, i] ** 2 + xy[1, i] ** 2) + np.sqrt(xy[1, i] ** 2 + (xy[0, i] - Ri) ** 2)) / 3e8
    return toa


if __name__ == '__main__':
    Rr = 4
    Ri = 20
    rho = np.pi / 20
    # ------------------theroratial curve-------------------
    c = 3e8
    stnum = 20
    tau_min = (2 * Rr + 2 * Ri * np.cos(rho) - 2 ** 0.5 * np.sqrt(
        -Ri ** 2 + 2 * Rr ** 2 + Ri ** 2 * np.cos(2 * rho))) / (2 * c)
    tau_max = (2 * Rr + 2 * Ri * np.cos(rho) + 2 ** 0.5 * np.sqrt(
        -Ri ** 2 + 2 * Rr ** 2 + Ri ** 2 * np.cos(2 * rho))) / (2 * c)
    tau = np.zeros(stnum)
    for i in range(0, stnum):
        tau[i] = tau_min + i * (tau_max - tau_min) / stnum
    A = (Rr ** 2) * (2 * np.arccos(np.sin(rho) * Ri / Rr) - np.sin(2 * np.arccos(np.sin(rho) * Ri / Rr)))
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
    # -------------------Actual curve-----------------------
    num = 50000
    xy = np.zeros((2, num))
    for i in range(0, num):
        xy[:, i] = UniCutXY(Ri, 0, Rr, rho, 1)
    toa = GetTOA(xy, Ri)
    # -------------------Draw curve-----------------------
    plt.plot(tau, f, color='k')
    plt.hist(toa, tau, density='true', color='#B0C4DE')
    plt.legend(["TOA pdf", "Normalized histogram"])
    plt.title('PDF of TOA')
    plt.xlabel('Time of Arrival(s)')
    plt.ylabel('Probability Density')
    plt.grid(True)
    plt.savefig("CutDiskTOAPDF.png")
    plt.show()
