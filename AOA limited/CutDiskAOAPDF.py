import numpy as np
import matplotlib.pyplot as plt


def UniCutXY(Rx, Ry, RadiusSqure, rho):
    while (1):
        theta = np.random.random() * 2 * np.pi
        R = np.random.uniform(0, RadiusSqure)
        x = np.sin(theta) * (R ** 0.5)
        y = np.cos(theta) * (R ** 0.5)
        # if (((y + Ry) >= ((x + Rx) * np.tan(rho))) | ((y + Ry) <= ((x + Rx) * np.tan(-rho)))):
        if (((y + Ry) >= ((x + Rx) * np.tan(rho)))):
            u = np.array((x + Rx, y + Ry))
            break
    return u


def GetAOA(xy):
    l = np.size(xy, 1)
    aoa = np.zeros(l)
    for i in range(0, l):
        aoa[i] = np.arctan(xy[1, i] / xy[0, i])
    return aoa


def GetTOA(xy, Ri):
    l = np.size(xy, 1)
    toa = np.zeros(l)
    for i in range(0, l):
        toa[i] = (np.sqrt(xy[0, i] ** 2 + xy[1, i] ** 2) + np.sqrt(xy[1, i] ** 2 + (xy[0, i] - Ri) ** 2)) / 2e8
    return toa


if __name__ == '__main__':
    Rr = 4
    Ri = 20
    rho = np.pi / 20
    # ------------------theroratial curve-------------------
    stnum = 20
    theta_max = np.arcsin(Rr / Ri)
    theta_b = np.zeros(stnum)
    for i in range(0, stnum):
        theta_b[i] = rho + i * (theta_max - rho) / stnum
    A = (Rr ** 2) * (2 * np.arccos(np.sin(rho) * Ri / Rr) - np.sin(2 * np.arccos(np.sin(rho) * Ri / Rr)))
    f = 2 * Ri * np.cos(theta_b) * np.sqrt((Ri ** 2) * (np.cos(theta_b) ** 2) - Ri ** 2 + Rr ** 2) / A
    f = f * 2
    # -------------------Actual curve-----------------------
    num = 50000
    xy = np.zeros((2, num))
    for i in range(0, num):
        xy[:, i] = UniCutXY(Ri, 0, Rr ** 2, rho)
    aoa = GetAOA(xy)
    toa = GetTOA(xy, Ri)
    # -------------------Draw curve-----------------------
    plt.plot(theta_b, f, color='k')
    plt.hist(aoa,theta_b, density='true', color='#B0C4DE')
    plt.legend(["AOA pdf", "Normalized histogram"])
    plt.title('PDF of AOA')
    plt.xlabel('Angel of Arrival(rad)')
    plt.ylabel('Probability Density')
    plt.grid(True)
    plt.savefig("CutDiskAOAPDF.png")
    plt.show()

