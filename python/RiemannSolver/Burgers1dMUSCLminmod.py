import numpy as np
#MUSCL + minmod for 1d Burgers equation
#介绍：https://zhuanlan.zhihu.com/p/331771977
nmax = 8
tmax = 100
dx = 1.0 / nmax
dt = 1.0 / tmax
U = np.zeros((tmax,nmax))
uL = np.zeros((tmax,nmax))#U_{i+1/2}^L
uR = np.zeros((tmax,nmax))#U_{i-1/2}^R
Fp = np.zeros((tmax,nmax))#F plus 1/2
Fm = np.zeros((tmax,nmax))#F minus 1/2
for i in range(0,nmax):
    if(i < nmax/2):
        U[0,i] = 2
    else:
        U[0,i] = 1

cfl = 0.8
dx = 0.1
dt = 0.01

def minmod(a,b):
    if(a * b > 0):
        if(abs(a) > abs(b)):
            return b
        else:
            return a
    else:
        return 0
        
def Flux(u):
    return 0.5 * u * u#Burgers 方程

for k in range(0,tmax - 1):
    for i in range(1,nmax-1):
        a = (U[k,i+1] - U[k,i])/dx
        b = (U[k,i] - U[k,i-1])/dx
        uL[k,i] = U[k,i] + 0.5 * dx * minmod(a,b)
        uR[k,i] = U[k,i] - 0.5 * dx * minmod(a,b)
    
    uL[k,0] = uL[k,1]
    uL[k,-1] = uL[k,-2]
    uR[k,0] = uR[k,1]
    uR[k,-1] = uR[k,-2]
    
    for i in range(1,nmax-1):
        Fp[k,i] = 0.5*(Flux(uL[k,i+1]) + Flux(uR[k,i])) - 0.5 * dx / dt * (uR[k,i+1] - uL[k,i])
        Fm[k,i] = 0.5*(Flux(uL[k,i]) + Flux(uR[k,i-1])) - 0.5 * dx / dt * (uR[k,i] - uL[k,i-1])
        U[k+1,i] = U[k,i] - dt/dx*(Fp[k,i] - Fm[k,i])
    
    U[k+1,0] = U[k+1,1]
    U[k+1,-1] = U[k+1,-2]