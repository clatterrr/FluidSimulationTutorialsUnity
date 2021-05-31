import numpy as np
import random
"""
D:\FluidSim\FluidSim\LevelSet\LevelSet-m2aster
Conservation
"""

Lx = 4
Ly = 4
nmax = 100
mmax = 100
dx = Lx / nmax
dy = Ly / mmax
dt = 0.01
eps = dx / 2
lamb = dt / dx
x = np.zeros((nmax))
y = np.zeros((mmax))
for i in range(0,nmax):
    x[i] = dx * i - Lx / 2
for j in range(0,mmax):
    y[j] = dy * j - Ly / 2
    
phi = np.zeros((nmax,mmax))
phis = np.zeros((nmax,mmax))
U = np.zeros((nmax,mmax))
V = np.zeros((nmax,mmax))
grad = np.zeros((nmax,mmax))
xNormal = np.zeros((nmax,mmax))
yNormal = np.zeros((nmax,mmax))
Fleft = np.zeros((nmax,mmax))
Fright = np.zeros((nmax,mmax))
Gleft = np.zeros((nmax,mmax))
Gright = np.zeros((nmax,mmax))

F = np.zeros((nmax,mmax))
sxleft = np.zeros((nmax,mmax))
sxright = np.zeros((nmax,mmax))
sybot = np.zeros((nmax,mmax))
sytop = np.zeros((nmax,mmax))

sx = np.zeros((nmax,mmax))
sy = np.zeros((nmax,mmax))

FrightFlux = np.zeros((nmax,mmax))
FleftFlux = np.zeros((nmax,mmax))
GtopFlux = np.zeros((nmax,mmax))
GbotFlux = np.zeros((nmax,mmax))

phiLeftPlus = np.zeros((nmax,mmax))
phiLeftMinus = np.zeros((nmax,mmax))
phiRightPlus = np.zeros((nmax,mmax))
phiRightMinus = np.zeros((nmax,mmax))
phiTopPlus = np.zeros((nmax,mmax))
phiTopMinus = np.zeros((nmax,mmax))
phiBotPlus = np.zeros((nmax,mmax))
phiBotMinus = np.zeros((nmax,mmax))


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
        for j in range(1,mmax-1):
            if i == 1:
                sxleft[i,j] = superbee((phi[i,j] - phi[i-1,j])/dx,(phi[i-1,j] - 0)/dx)
            else:
                sxleft[i,j] = superbee((phi[i,j] - phi[i-1,j])/dx,(phi[i-1,j] - phi[i-2,j])/dx)
            if i == nmax-2:
                sxright[i,j] = superbee((0-phi[i+1,j])/dx,(phi[i+1,j]-phi[i,j])/dx)
            else:
                sxright[i,j] = superbee((phi[i+2,j]-phi[i+1,j])/dx,(phi[i+1,j]-phi[i,j])/dx)
            sx[i,j] = superbee((phi[i+1,j]-phi[i,j])/dx,(phi[i,j]-phi[i-1,j])/dx)
            
            if j == 1:
                sybot[i,j] = superbee((phi[i,j]-phi[i,j-1])/dy,(phi[i,j-1]-0)/dy)
            else:
                sybot[i,j] = superbee((phi[i,j]-phi[i,j-1])/dy,(phi[i,j-1] - phi[i,j-2])/dy)
            if j == mmax-2:
                sytop[i,j] = superbee((0-phi[i,j+1])/dy,(phi[i,j+1]-phi[i,j])/dy)
            else:
                sytop[i,j] = superbee((phi[i,j+2]-phi[i,j+1])/dy,(phi[i,j+1]-phi[i,j])/dy)
            sy[i,j] = superbee((phi[i,j+1]-phi[i,j])/dy,(phi[i,j]-phi[i,j-1])/dy)
            
            phiLeftPlus[i,j] = phi[i,j] - sx[i,j] * dx / 2
            phiLeftMinus[i,j] = phi[i-1,j] + sxleft[i,j] * dx / 2
            phiRightPlus[i,j] = phi[i+1,j] - sxright[i,j] * dx / 2
            phiRightMinus[i,j] = phi[i,j] + sx[i,j] * dx / 2
            
            phiTopPlus[i,j] = phi[i,j+1] - sytop[i,j] * dy / 2
            phiTopMinus[i,j] = phi[i,j] + sy[i,j] * dy / 2
            phiBotPlus[i,j] = phi[i,j] - sy[i,j] * dy / 2
            phiBotMinus[i,j] = phi[i,j-1] + sybot[i,j] * dy / 2
            
            FrightFlux[i,j] = max(U[i,j],0)*phiRightMinus[i,j] + min(U[i,j],0)*phiRightPlus[i,j]
            FleftFlux[i,j] = max(U[i-1,j],0)*phiLeftMinus[i,j] + min(U[i-1,j],0)*phiLeftPlus[i,j]

            GtopFlux[i,j] = max(V[i,j],0)*phiTopMinus[i,j] + min(V[i,j],0)*phiTopPlus[i,j]
            GbotFlux[i,j] = max(V[i,j-1],0)*phiBotMinus[i,j] + min(V[i,j-1],0)*phiBotPlus[i,j]
            
            F[i,j] = - (FrightFlux[i,j] - FleftFlux[i,j] + GtopFlux[i,j] - GbotFlux[i,j])/dx
            
            phis[i,j] = phi[i,j] + dt * F[i,j]
            
for i in range(0,nmax):
    for j in range(0,mmax):
        phi[i,j] = np.sqrt((x[i] - 0.5)**2 + (y[j]-0.75)**2) - 0.15
        
        U[i,j] = dx * j + dx / 2 - Lx / 2
        V[i,j] =  - (dy * i + dy / 2 - Ly / 2)
        
phi = 1 / (1 + np.exp(phi/eps))
for t in range(0,100):
    # ConserveLevelSetEvolve
    for i in range(1,nmax-1):
        for j in range(1,mmax-1):
            grad[i,j] = np.sqrt((phi[i+1,j] - phi[i-1,j])**2 + (phi[i,j+1] - phi[i,j-1])**2)
            xNormal[i,j] = (phi[i+1,j] - phi[i-1,j])/grad[i,j]
            yNormal[i,j] = (phi[i,j+1] - phi[i,j-1])/grad[i,j]
    # 重新初始化
    for k in range(300):
        if t % 5 != 4:
            break
        f = phi * (1 - phi) * xNormal
        g = phi * (1 - phi) * yNormal
        diff = 0
        for i in range(1,nmax-1):
            for j in range(1,mmax-1):
                Fright[i,j] = (f[i,j] + f[i+1,j])/2 - (phi[i+1,j] - phi[i,j])/2
                Fleft[i,j] = (f[i,j] + f[i-1,j])/2 - (phi[i,j] - phi[i-1,j])/2
            
                Gright[i,j] = (g[i,j] + g[i,j+1])/2 - (phi[i,j+1] - phi[i,j])/2
                Gleft[i,j] = (g[i,j] + g[i,j-1])/2 - (phi[i,j] - phi[i,j-1])/2
        for i in range(0,nmax):
            for j in range(0,mmax):
                dphi = lamb * (Fright[i,j] - Fleft[i,j] + Gright[i,j] - Gleft[i,j])
                phi[i,j] = phi[i,j] - dphi
                diff += dphi*dphi/nmax/mmax
        if diff < 1e-3:
            break
    oldphi = phi.copy()
    RungeKutta2()
    phi = phis.copy()
    RungeKutta2()
    phi = (oldphi + phis)/2
    
    