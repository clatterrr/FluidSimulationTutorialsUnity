import numpy as np
import random
"""
D:\FluidSim\FluidSim\LevelSet\LevelSet-m2aster
Conservation
"""

Lx = 1
nmax = 100
dx = Lx / nmax
dt = 0.01
eps = dx / 2
lamb = dt / dx
x = np.zeros((nmax))
for i in range(0,nmax):
    x[i] = dx * i - Lx / 2
    
phi = np.zeros((nmax))
phis = np.zeros((nmax))
U = np.zeros((nmax))
Fleft = np.zeros((nmax))
Fright = np.zeros((nmax))



def superbee(x,y):
    L = 0
    if x * y > 0:
        if (abs(y) > 2*abs(x)) |(abs(y) < abs(x)/2):
            L = 2 * np.sign(x) * min(abs(x),abs(y))
        else:
            L = np.sign(x) * min(abs(x),abs(y))
    return L

def RungeKutta2():
    for i in range(1,nmax-1):
            sxleft = 0
            sxright = 0
            if i == 1:
                sxleft = superbee((phi[i] - phi[i-1])/dx,(phi[i-1] - 0)/dx)
            else:
                sxleft = superbee((phi[i] - phi[i-1])/dx,(phi[i-1] - phi[i-2])/dx)
            if i == nmax-2:
                sxright = superbee((0-phi[i+1])/dx,(phi[i+1]-phi[i])/dx)
            else:
                sxright = superbee((phi[i+2]-phi[i+1])/dx,(phi[i+1]-phi[i])/dx)
            sx = superbee((phi[i+1]-phi[i])/dx,(phi[i]-phi[i-1])/dx)
            
            
            phiLeftPlus = phi[i] - sx * dx / 2
            phiLeftMinus = phi[i-1] + sxleft * dx / 2
            phiRightPlus = phi[i+1] - sxright * dx / 2
            phiRightMinus = phi[i] + sx * dx / 2
            
            
            FrightFlux= max(U[i],0)*phiRightMinus + min(U[i],0)*phiRightPlus
            FleftFlux = max(U[i-1],0)*phiLeftMinus + min(U[i-1],0)*phiLeftPlus

            
            phis[i] = phi[i] - dt * (FrightFlux - FleftFlux)/dx