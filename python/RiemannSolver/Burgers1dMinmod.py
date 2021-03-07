import numpy as np
#minmod limiter for 1d Burgers equation
#介绍：https://zhuanlan.zhihu.com/p/331771977
nmax = 8
tmax = 1000
dx = 1.0 / nmax
dt = 1.0 / tmax
U = np.zeros((tmax,nmax))
F = np.zeros((tmax,nmax))
upred = np.zeros((tmax,nmax))
for i in range(0,nmax):
    if(i < 2):
        U[0,i] = 2
    else:
        U[0,i] = 1
F[0,:] = 0.5 * U[0,:] * U[0,:]

cfl = 0.8
dx = 0.1
dt = 0.01

def minmod(a,b):
    if(b == 0):
        return 1
    elif(a * b <= 0):
        return 0
    phi = min(a/b,1)
    return phi

for k in range(0,tmax - 1):
    for i in range(1,nmax-1):
        a = ((U[k,i+1] - dt / dx * F[k,i+1]) - (U[k,i] - dt / dx * F[k,i]))/dx
        b = ((U[k,i] - dt / dx * F[k,i]) - (U[k,i-1] - dt / dx * F[k,i-1]))/dx
        upred[k,i] = U[k,i] + 0.5 * dx * a *  minmod(a,b)
    
    upred[k,0] = upred[k,1]
    upred[k,-1] = upred[k,-2]#边界情况
    F[k+1,:] = 0.5 * upred[k,:] * upred[k,:]
            
    for i in range(1,nmax-1):
        U[k+1,i] = U[k,i] - dt/dx*(F[k+1,i] - F[k+1,i-1])
    
    U[k+1,0] = U[k+1,1]
    U[k+1,-1] = U[k+1,-2]
    F[k+1,:] = 0.5 * U[k+1,:] * U[k+1,:]