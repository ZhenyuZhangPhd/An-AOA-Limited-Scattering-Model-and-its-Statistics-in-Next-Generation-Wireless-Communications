import numpy as np
import math
from numpy.linalg import linalg
# import tdoa
import matplotlib.pyplot as plt


def chan(Bs, Ms, snr):
    r0 = np.array([np.zeros((7, 1))])
    r = np.array([np.zeros((6, 1))])
    for j in range(len(Bs[:, 0])):
        r0[0, j] = np.sqrt(((Bs[j, 0] - Ms[0]) ** 2) + ((Bs[j, 1] - Ms[1]) ** 2))
    for i in range(len(Bs[:, 0]) - 1):
        r[0, i] = r0[0, i + 1] - r0[0, 0]

    r = awgn(r, snr)
    Q = np.eye(len(Bs) - 1)
    k1 = 0
    k = np.array([np.zeros((6, 1))])
    for i in range(len(Bs[:, 0])-1):
        # k here si differen from matlab for the k[0 , 0] = 0 with matlab's = 7500
        k[0,i]= Bs[i+1,0] ** 2 + Bs[i+1,1] ** 2
    Ga=np.zeros((6, 3))
    for i in range(len(Bs[:, 0])-1):
        Ga[i,0]= -Bs[i+1,0]
        Ga[i, 1] = -Bs[i + 1, 1]
        Ga[i, 2] = -r[0,i]
    h=np.array([np.zeros((6,1))])
    for i in range(len(Bs[:, 0]) - 1):
        h[0,i]= 0.5*(r[0,i]**2 -k[0,i]+k1)
    ####
    Q=np.mat(Q)
    h = np.mat(h)
    Ga = np.mat(Ga)
    Za0=(Ga.T*Q.I*Ga).I*Ga.T*Q.I*h.T
    B=np.eye(len(Bs)-1)
    B=np.mat(B)
    for i in range(len(Bs[:, 0]) - 1):
        B[i,i]=np.sqrt((Bs[i, 0]-Za0[0,0])**2 +(Bs[i, 1]-Za0[1,0])**2)
    FI=B*Q*B
    Za1 = (Ga.T*FI.I*Ga).I*Ga.T*FI.I*h.T
    CovZa=(Ga.T*FI.I*Ga).I
    sB=np.eye(3)
    for i in range(3):
        sB[i,i] = Za1[i,0]
    sFI = 4*sB*CovZa*sB
    sGa = np.mat([[1,0],[0,1],[1,1]])
    sh=np.mat([[Za1[0,0]**2],[Za1[1,0]**2],[Za1[2,0]**2]])
    Za2=(sGa.T*sFI.I*sGa).I*sGa.T*sFI.I*sh
    sZ = np.sqrt(abs(Za2))
    xy_es = sZ
    return xy_es


def awgn(x, snr):
    snr = 10 ** (snr / 10.0)
    xpower = np.sum(x ** 2) / len(x)
    npower = xpower / snr
    n = np.random.randn(len(x)) * np.sqrt(npower)
    n = n.reshape(len(x), 1)
    y = x + n
    return y


def drawthepic():
    plt.show()


if __name__ == '__main__':
    Bs = np.array([[0, 0], [np.sqrt(3), 0], [0.5 * np.sqrt(3), -1.5], [-0.5 * np.sqrt(3), 1.5], [-np.sqrt(3), 0],
                   [-0.5 * np.sqrt(3), -1.5], [0.5 * np.sqrt(3), 1.5]])
    Bs = Bs*50

    fig1 = plt.figure()

    ax1=fig1.add_subplot(111, aspect='equal')
    # draw the seven base station
    plt.scatter(Bs[:, 0], Bs[:, 1])

    rect=plt.Rectangle(
        (-10, -10),  # (x,y)矩形左下角
        20,  # width长
        20,  # height宽
        color='maroon',
        alpha=0.5
    )
    ax1.add_patch(rect)
    plt.legend(["Ms","Bs"])
    Ms = np.array([150, 200])
    snr = 100
    Ms_es = chan(Bs, Ms, snr)
    print(Ms_es)
    rmse = abs(Ms_es - Ms)
    plt.show()
    fig1.savefig('fig1',dpi=90)
    print("finished")

    # ax1.add_patch(rect)

