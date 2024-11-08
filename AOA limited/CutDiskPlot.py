import numpy as np
import math
from numpy.linalg import linalg
# import tdoa
import matplotlib.pyplot as plt

def UniDiskPlot(Rx,Ry,radius,num):
    for i in range(1,num):
        theta = np.random.random()*2*np.pi
        R = np.random.uniform(0,radius)
        x = np.sin(theta) * (R**0.5)
        y = np.cos(theta) * (R**0.5)
        # plt.plot(x+Rx, y+Ry, '*',color='b')
        plt.scatter(x+Rx, y+Ry,color='b')

def UniCutPlot(Rx,Ry,radius,num,rho):
    for i in range(1,num):
        theta = np.random.random()*2*np.pi
        R = np.random.uniform(0,radius)
        x = np.sin(theta) * (R**0.5)
        y = np.cos(theta) * (R**0.5)
        if(((y + Ry)>=((x + Rx)*np.tan(rho)))|((y + Ry)<=((x + Rx)*np.tan(-rho)))):
            u=plt.scatter(x + Rx, y + Ry,  color='k')
    return u



if __name__ == '__main__':
    plt.figure(figsize=(10, 5), dpi=125)
    print("start\n")
    Rix=10
    Rrx=4
    rho=np.pi/12
    blockx=3
    blocky=np.tan(rho)*blockx
    plt.axis('equal')
    # plt.plot([-1, 20], [0, 0],'k')
    # plt.plot([0, 0], [-1, 5], 'k')
    plt.plot([0, 20], [0, 20*np.tan(rho)], '--',color='k')
    plt.plot([0, 20], [0, 20 * np.tan(-rho)], '--', color='k')
    Bs=plt.scatter(0,0)
    Ms=plt.scatter(Rix,0)
    plt.plot((blockx,blockx),(blocky,-blocky), color='k')
    x = np.linspace(6, 14)
    y = np.sqrt(16 - (x - 10) ** 2)  # 绘制上半圆
    plt.plot(x, y,'--',color='k')
    y = - np.sqrt(16 - (x - 10)**2)  # 绘制下半圆
    plt.plot(x, y,'--',color='k')
    u=UniCutPlot(Rix,0,16,800,rho)
    plt.legend([Bs, Ms,u], ["Bs", "Ms","st"])
    print("random")
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Random Scatter')
    plt.grid(True)
    plt.savefig("CutDiskScatterer.png")
    plt.show()

    print("finished")
